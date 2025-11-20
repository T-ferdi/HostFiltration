"""
Data Preprocessing Module
Handles data loading, tokenization, and dataset preparation for the token classification model.
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from transformers import AutoTokenizer
from datasets import Dataset


def load_data(csv_path="../subsequences_dataset.csv"):
    """Load and split the dataset into train, validation, and test sets."""
    df = pd.read_csv(csv_path)
    df["label"] = df["label"].astype(float)
    
    train_df, temp_df = train_test_split(
        df, test_size=0.3, stratify=df['label'], random_state=42
    )
    val_df, test_df = train_test_split(
        temp_df, test_size=0.5, stratify=temp_df['label'], random_state=42
    )
    
    print("Label distribution:")
    print(df['label'].value_counts())
    print(f"Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")
    
    return train_df, val_df, test_df


def tokenize(batch, tokenizer):
    """Tokenize sequences in a batch."""
    return tokenizer(batch["sequence"], padding="max_length", truncation=True, max_length=512)


def broadcast_labels(batch):
    """Create token-level labels by repeating the sequence label."""
    seq_len = len(batch['input_ids'][0])
    token_labels = []
    for l in batch['label']:
        token_labels.append([l] * seq_len)
    batch['labels'] = token_labels
    return batch


def prepare_datasets(train_df, val_df, test_df, tokenizer):
    """Prepare datasets for training by tokenizing and adding labels."""
    # Convert to HuggingFace datasets
    train_ds = Dataset.from_pandas(train_df)
    val_ds = Dataset.from_pandas(val_df)
    test_ds = Dataset.from_pandas(test_df)
    
    # Tokenize
    print("Tokenizing datasets...")
    train_ds = train_ds.map(lambda x: tokenize(x, tokenizer), batched=True)
    val_ds = val_ds.map(lambda x: tokenize(x, tokenizer), batched=True)
    test_ds = test_ds.map(lambda x: tokenize(x, tokenizer), batched=True)
    
    # Add labels
    print("Adding labels...")
    train_ds = train_ds.map(broadcast_labels, batched=True)
    val_ds = val_ds.map(broadcast_labels, batched=True)
    test_ds = test_ds.map(broadcast_labels, batched=True)
    
    # Set format for PyTorch
    train_ds.set_format(type="torch", columns=["input_ids", "attention_mask", "labels", "contig"])
    val_ds.set_format(type="torch", columns=["input_ids", "attention_mask", "labels", "contig"])
    test_ds.set_format(type="torch", columns=["input_ids", "attention_mask", "labels", "contig"])
    
    return train_ds, val_ds, test_ds