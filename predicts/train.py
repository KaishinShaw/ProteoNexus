#!/usr/bin/env python3
import sys
import ast
import pandas as pd
import numpy as np
import os
import pickle
import joblib
import json
import gc
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')
import traceback

if len(sys.argv) < 2:
    print("Usage: python3 my_script.py \"['COLUMN_NAME']\"")
    sys.exit(1)

target_columns = ast.literal_eval(sys.argv[1])
print("Column running:", target_columns)

proteo_df = pd.read_pickle('f_proteo_train.pkl')

# Cell 1: Import necessary libraries
# Machine Learning
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (roc_auc_score, roc_curve, accuracy_score, 
                           precision_score, recall_score, f1_score, 
                           confusion_matrix, classification_report)
import xgboost as xgb
import lightgbm as lgb

# Deep Learning
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

# Hyperparameter Optimization
import optuna
from optuna.visualization import plot_optimization_history, plot_param_importances

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import auc

# Set random seeds for reproducibility
np.random.seed(42)
torch.manual_seed(42)

# GPU Memory Management Function
def clear_gpu_memory():
    """Clear GPU memory and garbage collect"""
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.synchronize()

# Cell 2: Data Preparation Functions
def prepare_data(df, target_column):
    """
    Prepare data for training by separating features and target
    """
    # Identify PHD columns
    phd_columns = [col for col in df.columns if col.startswith('PHD')]
    
    # Features are all columns except PHD columns
    feature_columns = [col for col in df.columns if col not in phd_columns]
    
    # Extract features and target
    X = df[feature_columns].values
    y = df[target_column].values
    
    return X, y, feature_columns

def split_and_scale_data(X, y, test_size=0.2, random_state=42):
    """
    Split data into train/validation sets and apply standardization
    """
    # Check if we have enough samples for stratified split
    min_samples_per_class = 2
    unique_classes, class_counts = np.unique(y, return_counts=True)
    
    # Check if any class has too few samples
    if len(unique_classes) < 2:
        raise ValueError(f"Only one class present in the data. Need at least 2 classes for classification.")
    
    for cls, count in zip(unique_classes, class_counts):
        if count < min_samples_per_class:
            raise ValueError(f"Class {cls} has only {count} samples. Need at least {min_samples_per_class} samples per class.")
    
    # Check if we have enough samples for the split
    min_samples_for_split = int(1 / test_size) + 1
    if len(y) < min_samples_for_split:
        raise ValueError(f"Need at least {min_samples_for_split} samples to split with test_size={test_size}. Got {len(y)} samples.")
    
    # Split data
    X_train, X_val, y_train, y_val = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )
    
    # Standardize features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_val_scaled = scaler.transform(X_val)
    
    return X_train_scaled, X_val_scaled, y_train, y_val, scaler

# Cell 3: PyTorch MLP Model Definition
class MLPClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dims=[128, 64], dropout_rate=0.3):
        super(MLPClassifier, self).__init__()
        
        layers = []
        prev_dim = input_dim
        
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            ])
            prev_dim = hidden_dim
        
        layers.append(nn.Linear(prev_dim, 1))
        layers.append(nn.Sigmoid())
        
        self.model = nn.Sequential(*layers)
    
    def forward(self, x):
        return self.model(x)

def train_pytorch_model(model, train_loader, val_loader, epochs=100, lr=0.001, device='cuda'):
    """
    Train PyTorch model with early stopping
    """
    model = model.to(device)
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)
    
    best_val_loss = float('inf')
    patience = 10
    patience_counter = 0
    
    train_losses = []
    val_losses = []
    
    for epoch in range(epochs):
        # Training
        model.train()
        train_loss = 0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            
            optimizer.zero_grad()
            outputs = model(X_batch)
            outputs = outputs.squeeze(-1)
            loss = criterion(outputs, y_batch.float())
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                outputs = model(X_batch)
                outputs = outputs.squeeze(-1)
                loss = criterion(outputs, y_batch.float())
                val_loss += loss.item()
        
        avg_train_loss = train_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)
        
        train_losses.append(avg_train_loss)
        val_losses.append(avg_val_loss)
        
        # Early stopping
        if avg_val_loss < best_val_loss:
            best_val_loss = avg_val_loss
            patience_counter = 0
            best_model_state = model.state_dict()
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch}")
                break
    
    # Load best model
    model.load_state_dict(best_model_state)
    
    # Clear GPU memory after training
    del optimizer
    clear_gpu_memory()
    
    return model, train_losses, val_losses

# Cell 4: Hyperparameter Optimization Functions with Error Handling
def safe_optimize_function(optimize_func, trial, *args):
    """Wrapper function to safely execute optimization with error handling"""
    try:
        return optimize_func(trial, *args)
    except Exception as e:
        print(f"Error in optimization trial: {str(e)}")
        return 0.0  # Return worst possible score

def optimize_logistic_regression(trial, X_train, y_train, X_val, y_val):
    """Optuna optimization for Logistic Regression"""
    try:
        params = {
            'C': trial.suggest_float('C', 0.001, 10.0, log=True),
            'penalty': trial.suggest_categorical('penalty', ['l1', 'l2']),
            'solver': 'liblinear' if trial.params['penalty'] == 'l1' else 'lbfgs',
            'max_iter': 300,
            'n_jobs': -1
        }
        
        model = LogisticRegression(**params)
        model.fit(X_train, y_train)
        
        y_pred_proba = model.predict_proba(X_val)[:, 1]
        score = roc_auc_score(y_val, y_pred_proba)
        
        # Memory cleanup
        del model
        gc.collect()
        
        return score
    except Exception as e:
        print(f"Error in Logistic Regression optimization: {str(e)}")
        return 0.0

def optimize_xgboost(trial, X_train, y_train, X_val, y_val):
    """Optuna optimization for XGBoost with GPU"""
    try:
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 50, 300),
            'max_depth': trial.suggest_int('max_depth', 3, 9),
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
            'subsample': trial.suggest_float('subsample', 0.6, 1.0),
            'colsample_bytree': trial.suggest_float('colsample_bytree', 0.6, 1.0),
            'gamma': trial.suggest_float('gamma', 0.01, 1.0, log=True),
            'tree_method': 'gpu_hist' if torch.cuda.is_available() else 'hist',
            'predictor': 'cpu_predictor',
            'use_label_encoder': False,
            'eval_metric': 'logloss'
        }
        
        model = xgb.XGBClassifier(**params)
        model.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=False)
        
        y_pred_proba = model.predict_proba(X_val)[:, 1]
        score = roc_auc_score(y_val, y_pred_proba)
        
        # Memory cleanup
        del model
        clear_gpu_memory()
        
        return score
    except Exception as e:
        print(f"Error in XGBoost optimization: {str(e)}")
        return 0.0

def optimize_lightgbm(trial, X_train, y_train, X_val, y_val):
    """Optuna optimization for LightGBM with GPU"""
    try:
        params = {
            'n_estimators': trial.suggest_int('n_estimators', 50, 300),
            'num_leaves': trial.suggest_int('num_leaves', 20, 300),
            'max_depth': trial.suggest_int('max_depth', 3, 15),
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
            'feature_fraction': trial.suggest_float('feature_fraction', 0.5, 1.0),
            'bagging_fraction': trial.suggest_float('bagging_fraction', 0.5, 1.0),
        }
        
        if torch.cuda.is_available():
            params.update({
                'device': 'gpu',
                'gpu_use_dp': False,
            })
        
        model = lgb.LGBMClassifier(**params)
        model.fit(X_train, y_train, eval_set=[(X_val, y_val)], 
                  callbacks=[lgb.early_stopping(10), lgb.log_evaluation(0)])
        
        y_pred_proba = model.predict_proba(X_val)[:, 1]
        score = roc_auc_score(y_val, y_pred_proba)
        
        # Memory cleanup
        del model
        clear_gpu_memory()
        
        return score
    except Exception as e:
        print(f"Error in LightGBM optimization: {str(e)}")
        return 0.0

def optimize_mlp(trial, X_train, y_train, X_val, y_val, input_dim):
    """Optuna optimization for MLP"""
    try:
        # Hyperparameters
        n_layers = trial.suggest_int('n_layers', 1, 3)
        hidden_dims = []
        for i in range(n_layers):
            hidden_dims.append(trial.suggest_int(f'hidden_dim_{i}', 32, 256))
        
        dropout_rate = trial.suggest_float('dropout_rate', 0.1, 0.5)
        lr = trial.suggest_float('lr', 0.0001, 0.01, log=True)
        batch_size = trial.suggest_categorical('batch_size', [32, 64, 128])
        
        # Ensure batch_size is not larger than dataset size
        batch_size = min(batch_size, len(X_train))
        
        # Prepare data loaders
        train_dataset = TensorDataset(torch.FloatTensor(X_train), torch.LongTensor(y_train))
        val_dataset = TensorDataset(torch.FloatTensor(X_val), torch.LongTensor(y_val))
        
        # Set drop_last based on dataset size
        drop_last = len(train_dataset) > batch_size
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=drop_last)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        
        # Create and train model
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model = MLPClassifier(input_dim, hidden_dims, dropout_rate)
        
        trained_model, _, _ = train_pytorch_model(model, train_loader, val_loader, 
                                                 epochs=50, lr=lr, device=device)
        
        # Evaluate
        trained_model.eval()
        trained_model = trained_model.to('cpu')
        with torch.no_grad():
            y_pred_proba = trained_model(torch.FloatTensor(X_val)).squeeze(-1).numpy()
        
        score = roc_auc_score(y_val, y_pred_proba)
        
        # Memory cleanup
        del trained_model, model
        clear_gpu_memory()
        
        return score
    except Exception as e:
        print(f"Error in MLP optimization: {str(e)}")
        return 0.0

# Cell 5: Model Training and Evaluation Functions
def evaluate_model(model, X_val, y_val, model_name, is_pytorch=False):
    """
    Evaluate model performance and return metrics
    """
    try:
        if is_pytorch:
            model.eval()
            with torch.no_grad():
                y_pred_proba = model(torch.FloatTensor(X_val)).squeeze(-1).numpy()
        else:
            y_pred_proba = model.predict_proba(X_val)[:, 1]
        
        y_pred = (y_pred_proba >= 0.5).astype(int)
        
        metrics = {
            'model': model_name,
            'auc_roc': roc_auc_score(y_val, y_pred_proba),
            'accuracy': accuracy_score(y_val, y_pred),
            'precision': precision_score(y_val, y_pred, zero_division=0),
            'recall': recall_score(y_val, y_pred, zero_division=0),
            'f1': f1_score(y_val, y_pred, zero_division=0)
        }
        
        fpr, tpr, _ = roc_curve(y_val, y_pred_proba)
        
        return metrics, fpr, tpr, y_pred_proba
    except Exception as e:
        print(f"Error evaluating {model_name}: {str(e)}")
        # Return default metrics in case of error
        return {
            'model': model_name,
            'auc_roc': 0.5,
            'accuracy': 0.0,
            'precision': 0.0,
            'recall': 0.0,
            'f1': 0.0
        }, np.array([0, 1]), np.array([0, 1]), np.array([])

def plot_roc_curves(results_dict, target_column, save_path):
    """
    Plot ROC curves for all models with error handling
    """
    plt.figure(figsize=(10, 8))
    
    has_valid_models = False
    
    for model_type in ['optimized', 'default']:
        for result in results_dict.get(model_type, []):
            if 'fpr' in result and 'tpr' in result and len(result['fpr']) > 0:
                label = f"{result['model']} ({model_type})"
                plt.plot(result['fpr'], result['tpr'], 
                        label=f"{label} (AUC = {result['metrics']['auc_roc']:.3f})")
                has_valid_models = True
    
    if has_valid_models:
        plt.plot([0, 1], [0, 1], 'k--', label='Random')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title(f'ROC Curves - {target_column}')
        plt.legend(loc="lower right")
        plt.grid(True, alpha=0.3)
    else:
        plt.text(0.5, 0.5, 'No valid models to plot', 
                horizontalalignment='center', verticalalignment='center',
                transform=plt.gca().transAxes, fontsize=12)
        plt.title(f'ROC Curves - {target_column} (No valid models)')
    
    plt.savefig(os.path.join(save_path, f'roc_curves_{target_column}.png'), 
                dpi=300, bbox_inches='tight')
    plt.close()
    plt.clf()
    clear_gpu_memory()

def save_metrics_comparison(results_dict, target_column, save_path):
    """
    Save and visualize metrics comparison as separate plots with error handling
    """
    # Create comparison DataFrame
    all_metrics = []
    for model_type in ['default', 'optimized']:
        for result in results_dict.get(model_type, []):
            if 'metrics' in result:
                metrics = result['metrics'].copy()
                metrics['type'] = model_type
                all_metrics.append(metrics)
    
    if not all_metrics:
        print(f"No valid metrics to save for {target_column}")
        return pd.DataFrame()
    
    df_metrics = pd.DataFrame(all_metrics)
    
    # Save to CSV
    df_metrics.to_csv(os.path.join(save_path, f'metrics_comparison_{target_column}.csv'), 
                      index=False)
    
    # Create separate visualizations for each metric
    metrics_to_plot = ['auc_roc', 'accuracy', 'precision', 'f1']
    
    for metric in metrics_to_plot:
        plt.figure(figsize=(10, 6))
        
        try:
            df_pivot = df_metrics.pivot(index='model', columns='type', values=metric)
            if not df_pivot.empty:
                df_pivot.plot(kind='bar')
                plt.title(f'{metric.upper()} Comparison - {target_column}')
                plt.ylabel(metric)
                plt.xlabel('Model')
                plt.xticks(rotation=45, ha='right')
                plt.legend(title='Model Type')
                plt.grid(True, alpha=0.3)
            else:
                plt.text(0.5, 0.5, f'No valid data for {metric}', 
                        horizontalalignment='center', verticalalignment='center',
                        transform=plt.gca().transAxes, fontsize=12)
                plt.title(f'{metric.upper()} Comparison - {target_column} (No data)')
        except Exception as e:
            print(f"Error plotting {metric} for {target_column}: {str(e)}")
            plt.text(0.5, 0.5, f'Error plotting {metric}', 
                    horizontalalignment='center', verticalalignment='center',
                    transform=plt.gca().transAxes, fontsize=12)
            plt.title(f'{metric.upper()} Comparison - {target_column} (Error)')
        
        plt.tight_layout()
        
        # Save individual plot
        plt.savefig(os.path.join(save_path, f'{metric}_{target_column}.png'), 
                    dpi=300, bbox_inches='tight')
        plt.close()
        plt.clf()
    
    clear_gpu_memory()
    
    return df_metrics

# Cell 6: Independent Model Training Functions with Error Handling
def train_logistic_regression(X_train, y_train, X_val, y_val, target_column, 
                             optimize=True, n_trials=50, save_path='./models'):
    """
    Train Logistic Regression independently with error handling
    """
    try:
        model_save_path = os.path.join(save_path, target_column)
        os.makedirs(model_save_path, exist_ok=True)
        
        if optimize:
            # Optimize hyperparameters
            study = optuna.create_study(direction='maximize')
            study.optimize(lambda trial: safe_optimize_function(optimize_logistic_regression, trial, X_train, y_train, X_val, y_val), 
                          n_trials=n_trials, n_jobs=-1)
            
            # Train with best parameters
            best_params = study.best_params.copy()
            if best_params['penalty'] == 'l1':
                best_params['solver'] = 'liblinear'
            else:
                best_params['solver'] = 'lbfgs'
            best_params['max_iter'] = 300
            best_params['n_jobs'] = -1
            
            model = LogisticRegression(**best_params)
        else:
            # Use default parameters
            model = LogisticRegression(max_iter=300, n_jobs=-1)
        
        model.fit(X_train, y_train)
        
        # Evaluate
        metrics, fpr, tpr, y_pred_proba = evaluate_model(model, X_val, y_val, 'Logistic Regression')
        
        # Save model
        joblib.dump(model, os.path.join(model_save_path, 'logistic_regression_model.pkl'))
        
        print(f"Logistic Regression - AUC-ROC: {metrics['auc_roc']:.4f}")
        
        clear_gpu_memory()
        
        return model, metrics
        
    except Exception as e:
        print(f"Error training Logistic Regression for {target_column}: {str(e)}")
        traceback.print_exc()
        return None, None

def train_xgboost(X_train, y_train, X_val, y_val, target_column, 
                  optimize=True, n_trials=50, save_path='./models'):
    """
    Train XGBoost independently with error handling
    """
    try:
        model_save_path = os.path.join(save_path, target_column)
        os.makedirs(model_save_path, exist_ok=True)
        
        if optimize:
            # Optimize hyperparameters
            study = optuna.create_study(direction='maximize')
            study.optimize(lambda trial: safe_optimize_function(optimize_xgboost, trial, X_train, y_train, X_val, y_val), 
                          n_trials=n_trials)
            
            # Train with best parameters
            best_params = study.best_params.copy()
            best_params.update({
                'tree_method': 'gpu_hist' if torch.cuda.is_available() else 'hist',
                'predictor': 'cpu_predictor',
                'use_label_encoder': False,
                'eval_metric': 'logloss',
                'random_state': 42
            })
            
            model = xgb.XGBClassifier(**best_params)
        else:
            # Use default parameters
            model = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
        
        model.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=False)
        
        # Evaluate
        metrics, fpr, tpr, y_pred_proba = evaluate_model(model, X_val, y_val, 'XGBoost')
        
        # Save model
        joblib.dump(model, os.path.join(model_save_path, 'xgboost_model.pkl'))
        
        print(f"XGBoost - AUC-ROC: {metrics['auc_roc']:.4f}")
        
        clear_gpu_memory()
        
        return model, metrics
        
    except Exception as e:
        print(f"Error training XGBoost for {target_column}: {str(e)}")
        traceback.print_exc()
        return None, None

def train_lightgbm(X_train, y_train, X_val, y_val, target_column, 
                   optimize=True, n_trials=50, save_path='./models'):
    """
    Train LightGBM independently with error handling
    """
    try:
        model_save_path = os.path.join(save_path, target_column)
        os.makedirs(model_save_path, exist_ok=True)
        
        if optimize:
            # Optimize hyperparameters
            study = optuna.create_study(direction='maximize')
            study.optimize(lambda trial: safe_optimize_function(optimize_lightgbm, trial, X_train, y_train, X_val, y_val), 
                          n_trials=n_trials)
            
            # Train with best parameters
            best_params = study.best_params.copy()
            best_params['random_state'] = 42
            if torch.cuda.is_available():
                best_params.update({
                    'device': 'gpu',
                    'gpu_use_dp': False,
                })
            
            model = lgb.LGBMClassifier(**best_params)
        else:
            # Use default parameters
            model = lgb.LGBMClassifier(random_state=42)
        
        model.fit(X_train, y_train, eval_set=[(X_val, y_val)], 
                  callbacks=[lgb.early_stopping(10), lgb.log_evaluation(0)])
        
        # Set to CPU for inference if it was using GPU
        if 'device' in model.get_params() and model.get_params()['device'] == 'gpu':
            model.set_params(device='cpu')
        
        # Evaluate
        metrics, fpr, tpr, y_pred_proba = evaluate_model(model, X_val, y_val, 'LightGBM')
        
        # Save model
        joblib.dump(model, os.path.join(model_save_path, 'lightgbm_model.pkl'))
        
        print(f"LightGBM - AUC-ROC: {metrics['auc_roc']:.4f}")
        
        clear_gpu_memory()
        
        return model, metrics
        
    except Exception as e:
        print(f"Error training LightGBM for {target_column}: {str(e)}")
        traceback.print_exc()
        return None, None

def train_mlp(X_train, y_train, X_val, y_val, target_column, 
              optimize=True, n_trials=25, save_path='./models'):
    """
    Train MLP independently with error handling
    """
    try:
        model_save_path = os.path.join(save_path, target_column)
        os.makedirs(model_save_path, exist_ok=True)
        
        input_dim = X_train.shape[1]
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        if optimize:
            # Optimize hyperparameters
            study = optuna.create_study(direction='maximize')
            study.optimize(lambda trial: safe_optimize_function(optimize_mlp, trial, X_train, y_train, X_val, y_val, input_dim), 
                          n_trials=n_trials, n_jobs=-1)
            
            # Train with best parameters
            best_params = study.best_params
            n_layers = best_params['n_layers']
            hidden_dims = [best_params[f'hidden_dim_{i}'] for i in range(n_layers)]
            dropout_rate = best_params['dropout_rate']
            lr = best_params['lr']
            batch_size = best_params['batch_size']
        else:
            # Use default parameters
            hidden_dims = [128, 64]
            dropout_rate = 0.3
            lr = 0.001
            batch_size = 64
        
        # Ensure batch_size is not larger than dataset size
        batch_size = min(batch_size, len(X_train))
        
        # Prepare data loaders
        train_dataset = TensorDataset(torch.FloatTensor(X_train), torch.LongTensor(y_train))
        val_dataset = TensorDataset(torch.FloatTensor(X_val), torch.LongTensor(y_val))
        
        # Set drop_last based on dataset size
        drop_last = len(train_dataset) > batch_size
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, drop_last=drop_last)
        val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
        
        # Create and train model
        model = MLPClassifier(input_dim, hidden_dims, dropout_rate)
        model, train_losses, val_losses = train_pytorch_model(model, train_loader, val_loader, 
                                                              epochs=100, lr=lr, device=device)
        model = model.to('cpu')
        
        # Evaluate
        metrics, fpr, tpr, y_pred_proba = evaluate_model(model, X_val, y_val, 'MLP', is_pytorch=True)
        
        # Save model
        torch.save({
            'model_state_dict': model.state_dict(),
            'model_config': {
                'input_dim': input_dim,
                'hidden_dims': hidden_dims,
                'dropout_rate': dropout_rate
            }
        }, os.path.join(model_save_path, 'mlp_model.pth'))
        
        print(f"MLP - AUC-ROC: {metrics['auc_roc']:.4f}")
        
        clear_gpu_memory()
        
        return model, metrics
        
    except Exception as e:
        print(f"Error training MLP for {target_column}: {str(e)}")
        traceback.print_exc()
        return None, None

# Cell 7: Main Training Pipeline with Comprehensive Error Handling
def train_all_models(X_train, y_train, X_val, y_val, target_column, n_trials=50):
    """
    Train all models with default and optimized hyperparameters with error handling
    """
    results = {'default': [], 'optimized': []}
    trained_models = {}
    
    print(f"\nTraining models for {target_column}...")
    print("="*50)
    
    # Track failed models
    failed_models = []
    
    # 1. Logistic Regression
    print("\n1. Logistic Regression")
    try:
        # Default
        lr_default = LogisticRegression(max_iter=300, n_jobs=-1)
        lr_default.fit(X_train, y_train)
        metrics, fpr, tpr, _ = evaluate_model(lr_default, X_val, y_val, 'Logistic Regression')
        results['default'].append({'model': 'Logistic Regression', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
        print(f"Default AUC-ROC: {metrics['auc_roc']:.4f}")
        
        # Optimized
        lr_optimized, metrics_opt = train_logistic_regression(X_train, y_train, X_val, y_val, 
                                                             target_column, optimize=False, n_trials=n_trials)
        if lr_optimized is not None:
            metrics, fpr, tpr, _ = evaluate_model(lr_optimized, X_val, y_val, 'Logistic Regression')
            results['optimized'].append({'model': 'Logistic Regression', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
            trained_models['logistic_regression'] = lr_optimized
        else:
            failed_models.append('Logistic Regression (optimized)')
        
        # Memory cleanup
        del lr_default
        clear_gpu_memory()
        
    except Exception as e:
        print(f"Failed to train Logistic Regression: {str(e)}")
        failed_models.append('Logistic Regression')
    
    # 2. XGBoost
    print("\n2. XGBoost")
    try:
        # Default
        xgb_default = xgb.XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
        xgb_default.fit(X_train, y_train)
        metrics, fpr, tpr, _ = evaluate_model(xgb_default, X_val, y_val, 'XGBoost')
        results['default'].append({'model': 'XGBoost', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
        print(f"Default AUC-ROC: {metrics['auc_roc']:.4f}")
        
        # Optimized
        xgb_optimized, metrics_opt = train_xgboost(X_train, y_train, X_val, y_val, 
                                                  target_column, optimize=True, n_trials=n_trials)
        if xgb_optimized is not None:
            metrics, fpr, tpr, _ = evaluate_model(xgb_optimized, X_val, y_val, 'XGBoost')
            results['optimized'].append({'model': 'XGBoost', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
            trained_models['xgboost'] = xgb_optimized
        else:
            failed_models.append('XGBoost (optimized)')
        
        # Memory cleanup
        del xgb_default
        clear_gpu_memory()
        
    except Exception as e:
        print(f"Failed to train XGBoost: {str(e)}")
        failed_models.append('XGBoost')
    
    # 3. LightGBM
    print("\n3. LightGBM")
    try:
        # Default
        lgb_default = lgb.LGBMClassifier(random_state=42)
        lgb_default.fit(X_train, y_train)
        metrics, fpr, tpr, _ = evaluate_model(lgb_default, X_val, y_val, 'LightGBM')
        results['default'].append({'model': 'LightGBM', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
        print(f"Default AUC-ROC: {metrics['auc_roc']:.4f}")
        
        # Optimized
        lgb_optimized, metrics_opt = train_lightgbm(X_train, y_train, X_val, y_val, 
                                                   target_column, optimize=True, n_trials=n_trials)
        if lgb_optimized is not None:
            metrics, fpr, tpr, _ = evaluate_model(lgb_optimized, X_val, y_val, 'LightGBM')
            results['optimized'].append({'model': 'LightGBM', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
            trained_models['lightgbm'] = lgb_optimized
        else:
            failed_models.append('LightGBM (optimized)')
        
        # Memory cleanup
        del lgb_default
        clear_gpu_memory()
        
    except Exception as e:
        print(f"Failed to train LightGBM: {str(e)}")
        failed_models.append('LightGBM')
    
    # 4. MLP
    print("\n4. MLP (PyTorch)")
    try:
        # Default
        input_dim = X_train.shape[1]
        
        # Ensure batch size is appropriate for dataset size
        default_batch_size = min(64, len(X_train) // 2)
        
        train_dataset = TensorDataset(torch.FloatTensor(X_train), torch.LongTensor(y_train))
        val_dataset = TensorDataset(torch.FloatTensor(X_val), torch.LongTensor(y_val))
        
        # Set drop_last based on dataset size
        drop_last = len(train_dataset) > default_batch_size
        train_loader = DataLoader(train_dataset, batch_size=default_batch_size, shuffle=True, drop_last=drop_last)
        val_loader = DataLoader(val_dataset, batch_size=default_batch_size, shuffle=False)
        
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        print(f"Training MLP on {device}")
        
        mlp_default = MLPClassifier(input_dim, hidden_dims=[128, 64])
        mlp_default, _, _ = train_pytorch_model(mlp_default, train_loader, val_loader, device=device)
        mlp_default = mlp_default.to('cpu')
        
        metrics, fpr, tpr, _ = evaluate_model(mlp_default, X_val, y_val, 'MLP', is_pytorch=True)
        results['default'].append({'model': 'MLP', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
        print(f"Default AUC-ROC: {metrics['auc_roc']:.4f}")
        
        # Optimized
        mlp_optimized, metrics_opt = train_mlp(X_train, y_train, X_val, y_val, 
                                              target_column, optimize=True, n_trials=n_trials//2)
        if mlp_optimized is not None:
            metrics, fpr, tpr, _ = evaluate_model(mlp_optimized, X_val, y_val, 'MLP', is_pytorch=True)
            results['optimized'].append({'model': 'MLP', 'metrics': metrics, 'fpr': fpr, 'tpr': tpr})
            trained_models['mlp'] = mlp_optimized
        else:
            failed_models.append('MLP (optimized)')
        
        # Memory cleanup
        del mlp_default
        clear_gpu_memory()
        
    except Exception as e:
        print(f"Failed to train MLP: {str(e)}")
        failed_models.append('MLP')
    
    # Report failed models
    if failed_models:
        print(f"\nWarning: {len(failed_models)} models failed to train: {', '.join(failed_models)}")
    
    return results, trained_models

# Cell 8: Model Saving and Loading Functions with Error Handling
def save_models_and_results(trained_models, scaler, feature_columns, 
                            target_column, results, base_path='./models'):
    """
    Save trained models, scaler, and results with error handling
    """
    try:
        # Create directory for target column
        target_path = os.path.join(base_path, target_column)
        os.makedirs(target_path, exist_ok=True)
        
        # Save scaler
        joblib.dump(scaler, os.path.join(target_path, 'scaler.pkl'))
        
        # Save feature columns
        with open(os.path.join(target_path, 'feature_columns.json'), 'w') as f:
            json.dump(feature_columns, f)
        
        # Save results
        with open(os.path.join(target_path, 'training_results.json'), 'w') as f:
            # Convert numpy arrays to lists for JSON serialization
            results_json = {}
            for key in results:
                results_json[key] = []
                for item in results[key]:
                    item_copy = item.copy()
                    if 'fpr' in item_copy and hasattr(item_copy['fpr'], 'tolist'):
                        item_copy['fpr'] = item_copy['fpr'].tolist()
                    if 'tpr' in item_copy and hasattr(item_copy['tpr'], 'tolist'):
                        item_copy['tpr'] = item_copy['tpr'].tolist()
                    results_json[key].append(item_copy)
            json.dump(results_json, f)
        
        print(f"\nModels and results saved to {target_path}")
        
    except Exception as e:
        print(f"Error saving models and results for {target_column}: {str(e)}")
        traceback.print_exc()

def load_model_for_inference(target_column, model_type='best', base_path='./models'):
    """
    Load saved model for inference with error handling
    """
    try:
        target_path = os.path.join(base_path, target_column)
        
        # Load scaler and feature columns
        scaler = joblib.load(os.path.join(target_path, 'scaler.pkl'))
        with open(os.path.join(target_path, 'feature_columns.json'), 'r') as f:
            feature_columns = json.load(f)
        
        # Load training results to find best model
        with open(os.path.join(target_path, 'training_results.json'), 'r') as f:
            results = json.load(f)
        
        if model_type == 'best':
            # Find model with highest AUC-ROC
            best_auc = 0
            best_model_name = None
            for result in results.get('optimized', []):
                if result['metrics']['auc_roc'] > best_auc:
                    best_auc = result['metrics']['auc_roc']
                    best_model_name = result['model'].lower().replace(' ', '_')
            
            if best_model_name is None:
                raise ValueError(f"No valid models found for {target_column}")
            
            model_type = best_model_name
            print(f"Best model for {target_column}: {model_type} (AUC-ROC: {best_auc:.4f})")
        
        # Load the specified model
        if model_type == 'mlp':
            checkpoint = torch.load(os.path.join(target_path, f'{model_type}_model.pth'))
            model = MLPClassifier(**checkpoint['model_config'])
            model.load_state_dict(checkpoint['model_state_dict'])
            model.eval()
        else:
            model = joblib.load(os.path.join(target_path, f'{model_type}_model.pkl'))
        
        return model, scaler, feature_columns, model_type
        
    except Exception as e:
        print(f"Error loading model for {target_column}: {str(e)}")
        traceback.print_exc()
        return None, None, None, None

# Cell 9: Inference Function with Error Handling
def predict_probabilities(csv_path, target_column, model_type='best', base_path='./models'):
    """
    Load new data and predict probabilities with error handling
    """
    try:
        # Load model and preprocessing info
        model, scaler, feature_columns, selected_model = load_model_for_inference(
            target_column, model_type, base_path
        )
        
        if model is None:
            print(f"Failed to load model for {target_column}")
            return None
        
        # Load new data
        new_data = pd.read_csv(csv_path)
        
        # Ensure all required features are present
        missing_features = set(feature_columns) - set(new_data.columns)
        if missing_features:
            raise ValueError(f"Missing features in input data: {missing_features}")
        
        # Extract and order features correctly
        X_new = new_data[feature_columns].values
        
        # Scale features
        X_new_scaled = scaler.transform(X_new)
        
        # Make predictions
        if selected_model == 'mlp':
            with torch.no_grad():
                probabilities = model(torch.FloatTensor(X_new_scaled)).squeeze(-1).numpy()
        else:
            probabilities = model.predict_proba(X_new_scaled)[:, 1]
        
        # Create results DataFrame
        results_df = pd.DataFrame({
            'sample_id': range(len(probabilities)),
            'probability_positive': probabilities,
            'probability_negative': 1 - probabilities,
            'predicted_class': (probabilities >= 0.5).astype(int)
        })
        
        # Save results
        output_path = f'predictions_{target_column}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.csv'
        results_df.to_csv(output_path, index=False)
        
        print(f"\nPredictions saved to {output_path}")
        print(f"Model used: {selected_model}")
        print(f"Number of samples: {len(results_df)}")
        print(f"Positive predictions: {(results_df['predicted_class'] == 1).sum()}")
        print(f"Negative predictions: {(results_df['predicted_class'] == 0).sum()}")
        
        clear_gpu_memory()
        
        return results_df
        
    except Exception as e:
        print(f"Error making predictions for {target_column}: {str(e)}")
        traceback.print_exc()
        return None

# Main execution code
# Check if GPU is available
print(f"GPU Available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"GPU Device: {torch.cuda.get_device_name(0)}")

# Process each target column
all_results_file = 'all_results.pkl'
if os.path.exists(all_results_file):
    with open(all_results_file, 'rb') as f:
        all_results = pickle.load(f)
    print("Loaded existing all_results.")
else:
    all_results = {}

# Create summary for skipped columns
skipped_columns = []
failed_columns = []

for target_column in target_columns:
    
    if os.path.exists(all_results_file):
        with open(all_results_file, 'rb') as f:
            all_results = pickle.load(f)
    
    print(f"\n{'='*60}")
    print(f"Processing target: {target_column}")
    print(f"{'='*60}")
    
    try:
        # Drop NaN values for the target column
        proteo_dm = proteo_df.dropna(subset=[target_column])
        proteo = proteo_dm
        
        # Check if we have enough samples
        n_samples = len(proteo)
        if n_samples == 0:
            print(f"WARNING: No valid samples found for {target_column} (all values are NaN). Skipping this target.")
            skipped_columns.append({
                'target': target_column,
                'reason': 'No valid samples (all NaN)',
                'n_samples': 0
            })
            continue
        
        # Check class distribution
        unique_classes, class_counts = np.unique(proteo[target_column].values, return_counts=True)
        print(f"Total samples after dropping NaN: {n_samples}")
        print(f"Class distribution: {dict(zip(unique_classes, class_counts))}")
        
        # Check if we have both classes
        if len(unique_classes) < 2:
            print(f"WARNING: Only one class present for {target_column}. Need at least 2 classes for classification. Skipping this target.")
            skipped_columns.append({
                'target': target_column,
                'reason': f'Only one class present ({unique_classes[0]})',
                'n_samples': n_samples,
                'class_distribution': dict(zip(unique_classes, class_counts))
            })
            continue
        
        # Check minimum samples per class
        min_class_samples = min(class_counts)
        if min_class_samples < 2:
            print(f"WARNING: Insufficient samples in minority class for {target_column}. Skipping this target.")
            skipped_columns.append({
                'target': target_column,
                'reason': f'Insufficient samples in minority class (min={min_class_samples})',
                'n_samples': n_samples,
                'class_distribution': dict(zip(unique_classes, class_counts))
            })
            continue
        
        # Check if we have enough samples for train/test split
        min_samples_needed = 10
        if n_samples < min_samples_needed:
            print(f"WARNING: Insufficient total samples for {target_column} ({n_samples} < {min_samples_needed}). Skipping this target.")
            skipped_columns.append({
                'target': target_column,
                'reason': f'Insufficient total samples ({n_samples} < {min_samples_needed})',
                'n_samples': n_samples,
                'class_distribution': dict(zip(unique_classes, class_counts))
            })
            continue
        
        # Prepare data
        X, y, feature_columns = prepare_data(proteo, target_column)
        
        # Split and scale data
        X_train, X_val, y_train, y_val, scaler = split_and_scale_data(X, y)
        
        print(f"Training samples: {len(X_train)}")
        print(f"Validation samples: {len(X_val)}")
        print(f"Number of features: {len(feature_columns)}")
        print(f"Class distribution - Train: {np.bincount(y_train)}")
        print(f"Class distribution - Val: {np.bincount(y_val)}")
        
        # Train all models
        results, trained_models = train_all_models(X_train, y_train, X_val, y_val, 
                                                  target_column, n_trials=20)
        
        # Create save directory
        save_path = f'./results/{target_column}'
        os.makedirs(save_path, exist_ok=True)
        
        # Check if we have any valid results
        if any(results.values()):
            # Plot ROC curves
            plot_roc_curves(results, target_column, save_path)
            
            # Save metrics comparison
            metrics_df = save_metrics_comparison(results, target_column, save_path)
            
            # Save models and results
            save_models_and_results(trained_models, scaler, feature_columns, target_column, results)
            
            # Store results
            all_results[target_column] = {
                'results': results,
                'metrics_df': metrics_df,
                'trained_models': trained_models
            }
            
            with open(all_results_file, 'wb') as f:
                pickle.dump(all_results, f)
            print(f"Saved all_results for target '{target_column}' to '{all_results_file}'.")
        else:
            print(f"No valid models trained for {target_column}")
            failed_columns.append({
                'target': target_column,
                'reason': 'No models successfully trained'
            })
        
    except Exception as e:
        print(f"ERROR processing {target_column}: {str(e)}")
        traceback.print_exc()
        failed_columns.append({
            'target': target_column,
            'reason': f'Error during processing: {str(e)}',
            'n_samples': n_samples if 'n_samples' in locals() else 'Unknown'
        })
        continue
    
    # Clear memory after each target column
    clear_gpu_memory()
    gc.collect()

# Save summary of skipped columns
if skipped_columns:
    skipped_df = pd.DataFrame(skipped_columns)
    skipped_df.to_csv('skipped_columns_summary.csv', index=False)
    print(f"\n\nSkipped {len(skipped_columns)} columns. Summary saved to 'skipped_columns_summary.csv'")
    print("\nSkipped columns:")
    for skip_info in skipped_columns:
        print(f"  - {skip_info['target']}: {skip_info['reason']}")

# Save summary of failed columns
if failed_columns:
    failed_df = pd.DataFrame(failed_columns)
    failed_df.to_csv('failed_columns_summary.csv', index=False)
    print(f"\n\nFailed to process {len(failed_columns)} columns. Summary saved to 'failed_columns_summary.csv'")
    print("\nFailed columns:")
    for fail_info in failed_columns:
        print(f"  - {fail_info['target']}: {fail_info['reason']}")

print("\n\nTraining completed for all target columns!")

# Cell 11: Additional Analysis and Visualization Functions with Error Handling
def create_model_summary_report(all_results, save_path='./results'):
    """
    Create a comprehensive summary report for all models and targets with error handling
    """
    try:
        summary_data = []
        
        for target_column, target_results in all_results.items():
            try:
                if 'metrics_df' not in target_results or target_results['metrics_df'] is None:
                    print(f"No metrics found for {target_column}, skipping...")
                    continue
                
                metrics_df = target_results['metrics_df']
                
                if isinstance(metrics_df, pd.DataFrame) and not metrics_df.empty:
                    # Filter out any rows with NaN AUC-ROC values
                    valid_metrics = metrics_df.dropna(subset=['auc_roc'])
                    
                    if not valid_metrics.empty:
                        # Get best model for this target
                        best_model = valid_metrics.loc[valid_metrics['auc_roc'].idxmax()]
                        
                        # Calculate improvement if both default and optimized results exist
                        optimized_results = valid_metrics[valid_metrics['type'] == 'optimized']
                        default_results = valid_metrics[valid_metrics['type'] == 'default']
                        
                        improvement = 0.0
                        if not optimized_results.empty and not default_results.empty:
                            improvement = optimized_results['auc_roc'].max() - default_results['auc_roc'].max()
                        
                        summary_data.append({
                            'target': target_column,
                            'best_model': best_model['model'],
                            'best_auc_roc': best_model['auc_roc'],
                            'model_type': best_model['type'],
                            'improvement': improvement,
                            'n_models_trained': len(valid_metrics)
                        })
                    else:
                        print(f"No valid metrics for {target_column}")
                        
            except Exception as e:
                print(f"Error processing results for {target_column}: {str(e)}")
                continue
        
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_df = summary_df.sort_values('best_auc_roc', ascending=False)
            summary_df.to_csv(os.path.join(save_path, 'model_summary_report.csv'), index=False)
            
            print("\nModel Summary Report:")
            print(summary_df)
            
            # Create visualization of best models
            plt.figure(figsize=(12, 8))
            summary_df.plot(x='target', y='best_auc_roc', kind='bar')
            plt.title('Best AUC-ROC Scores by Target')
            plt.xlabel('Target Column')
            plt.ylabel('AUC-ROC Score')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(os.path.join(save_path, 'best_models_comparison.png'), dpi=300, bbox_inches='tight')
            plt.close()
            
            return summary_df
        else:
            print("No valid summary data to report")
            return pd.DataFrame()
            
    except Exception as e:
        print(f"Error creating model summary report: {str(e)}")
        traceback.print_exc()
        return pd.DataFrame()

# Generate summary report after training all models
if all_results:
    summary_report = create_model_summary_report(all_results)
else:
    print("\nNo results available to create summary report")