"""
Model Training Module
Handles model setup, training configuration, and training execution.
"""

import numpy as np
import torch
from sklearn.metrics import accuracy_score, f1_score
from transformers import (
    AutoConfig,
    AutoModelForTokenClassification,
    TrainingArguments,
    Trainer,
    TrainerCallback
)


def setup_device():
    """Setup and return the appropriate device (MPS, CUDA, or CPU)."""
    device = torch.device("mps" if torch.backends.mps.is_available() else 
                         "cuda" if torch.cuda.is_available() else "cpu")
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    print(f"MPS available: {torch.backends.mps.is_available()}")
    print(f"Using device: {device}")
    return device


def setup_model(device, model_name="InstaDeepAI/nucleotide-transformer-v2-50m-multi-species"):
    """Setup and configure the model for token classification."""
    config = AutoConfig.from_pretrained(model_name)
    config.num_labels = 2
    config.problem_type = "single_label_classification"
    
    model = AutoModelForTokenClassification.from_config(config)
    model.to(device)
    
    return model, config


def load_model_from_checkpoint(checkpoint_path, device):
    """Load a trained model from a checkpoint."""
    config = AutoConfig.from_pretrained(checkpoint_path)
    model = AutoModelForTokenClassification.from_config(config)
    model.to(device)
    return model, config


class LossLoggingCallback(TrainerCallback):
    """Callback to log training loss during training."""
    def on_log(self, args, state, control, **kwargs):
        if 'train_loss' in kwargs.get('logs', {}):
            print(f"Step {state.global_step}: Train Loss = {kwargs['logs']['train_loss']:.4f}")


def compute_metrics(p):
    """Compute accuracy and F1 score for sequence-level predictions."""
    # logits shape: (batch, seq_len, num_labels)
    logits = p.predictions
    probs = np.exp(logits) / np.exp(logits).sum(-1, keepdims=True)  # softmax

    # token predictions (0=bacteria, 1=human)
    token_preds = probs.argmax(-1)  # shape: (batch, seq_len)

    # Compute the percentage of human tokens per sequence
    human_token_pct = token_preds.sum(axis=1) / token_preds.shape[1]

    # Sequence-level predictions using 20% rule
    seq_preds = (human_token_pct >= 0.2).astype(int)

    # Sequence-level ground truth: first token label
    seq_labels = p.label_ids[:, 0].astype(int)

    # Metrics
    acc = accuracy_score(seq_labels, seq_preds)
    f1 = f1_score(seq_labels, seq_preds, average='binary')

    # Identify failures
    failures = []
    for i, (pred, label, pct, tokens) in enumerate(zip(seq_preds, seq_labels, human_token_pct, token_preds)):
        if pred != label:
            failures.append({
                "sequence_index": i,
                "predicted": int(pred),
                "true_label": int(label),
                "human_token_pct": float(pct),
                "token_preds": tokens.tolist()
            })

    return {
        "accuracy": acc,
        "f1": f1,
        "failures": failures
    }


def train_model(model, train_ds, val_ds, output_dir="./results", 
                learning_rate=1e-5, num_epochs=5, batch_size=8):
    """Train the model using HuggingFace Trainer."""
    training_args = TrainingArguments(
        output_dir=output_dir,
        eval_strategy="epoch",
        save_strategy="epoch",
        logging_dir="./logs",
        per_device_train_batch_size=batch_size,
        per_device_eval_batch_size=batch_size,
        num_train_epochs=num_epochs,
        weight_decay=0.01,
        load_best_model_at_end=True,
        learning_rate=learning_rate,
        warmup_steps=400,
        lr_scheduler_type="cosine",
        metric_for_best_model="f1",
        greater_is_better=True,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_ds,
        eval_dataset=val_ds,
        compute_metrics=compute_metrics,
        callbacks=[LossLoggingCallback()]
    )
    
    print("Starting training...")
    trainer.train()
    
    return trainer