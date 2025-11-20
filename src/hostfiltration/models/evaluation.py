"""
Model Evaluation Module
Handles model evaluation, testing, and visualization of results.
"""

import numpy as np
import torch
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


def get_sequence_labels(label_batch):
    """Safely extract sequence-level labels from token-level labels."""
    seq_labels = []
    for labels in label_batch:
        valid_labels = labels[labels != -100]  # ignore padding or special tokens
        if len(valid_labels) == 0:
            seq_labels.append(0)  # fallback
        else:
            seq_labels.append(int(valid_labels[0]))  # assume all valid tokens share same label
    return np.array(seq_labels)


def test_model_on_dataset(model, ds, batch_size=32, threshold=0.1):
    """Test the model on a dataset and return metrics."""
    model.eval()
    test_loader = torch.utils.data.DataLoader(ds, batch_size=batch_size, shuffle=False)
    all_preds = []
    all_labels = []
    all_attention_weights = None

    with torch.no_grad():
        for i, batch in enumerate(tqdm(test_loader, desc='Testing')):
            input_ids = batch['input_ids'].to(model.device)
            attention_mask = batch['attention_mask'].to(model.device)
            labels = batch['labels'].cpu().numpy()

            # Forward pass
            outputs = model(
                input_ids=input_ids,
                attention_mask=attention_mask,
                output_attentions=True
            )

            # Extract attention weights (if available)
            attention_weights = getattr(outputs, 'attentions', None)

            # Compute softmax probabilities
            logits = outputs.logits.cpu().numpy()
            probs = np.exp(logits) / np.exp(logits).sum(-1, keepdims=True)

            # Token-level predictions (0=bacteria, 1=human)
            token_preds = probs.argmax(-1)

            # Sequence-level predictions using threshold rule
            seq_preds = (token_preds.sum(axis=1) >= threshold * token_preds.shape[1]).astype(int)

            # True sequence labels
            seq_labels = get_sequence_labels(labels)

            all_preds.extend(seq_preds)
            all_labels.extend(seq_labels)

            # Save attention weights from first batch only
            if i == 0 and attention_weights is not None:
                all_attention_weights = [att.cpu().numpy() for att in attention_weights]

    # Compute metrics
    acc = accuracy_score(all_labels, all_preds)
    f1 = f1_score(all_labels, all_preds, average='binary')
    cm = confusion_matrix(all_labels, all_preds)

    print("\n--- Test Results ---")
    print(f"Test Accuracy: {acc:.4f}")
    print(f"Test F1 Score: {f1:.4f}")
    print("Confusion Matrix:")
    print(cm)

    # Identify failed cases
    failures = [
        {"index": i, "pred": p, "true": t}
        for i, (p, t) in enumerate(zip(all_preds, all_labels)) if p != t
    ]
    print(f"\nTotal failed sequences: {len(failures)} / {len(all_labels)}")

    return acc, f1, cm, all_attention_weights, failures


def evaluate_with_trainer(trainer, dataset):
    """Evaluate using the Trainer's predict method."""
    eval_results = trainer.predict(dataset)
    
    logits = eval_results.predictions
    labels = eval_results.label_ids

    # Softmax
    probs = np.exp(logits) / np.exp(logits).sum(-1, keepdims=True)

    # Token-level predictions
    token_preds = probs.argmax(-1)  # shape (batch, seq_len)

    # Sequence-level predictions (20% rule)
    seq_preds = (token_preds.sum(axis=1) >= 0.2 * token_preds.shape[1]).astype(int)

    # True labels (robust way to get them)
    seq_labels = np.array([lbl[lbl != -100][0] for lbl in labels])  # skip padding

    fail_indices = np.where(seq_preds != seq_labels)[0]
    success_indices = np.where(seq_preds == seq_labels)[0]

    print(f"Total sequences: {len(seq_labels)}")
    print(f"Failures: {len(fail_indices)}")
    print(f"Successes: {len(success_indices)}")
    
    return seq_preds, seq_labels, fail_indices, success_indices


def visualize_attention_heatmap(attention_weights, tokenizer, input_ids, layer_idx=-1, 
                                head_idx=0, sample_idx=0, max_tokens=50, figsize=(12, 10)):
    """Visualize attention weights as a heatmap."""
    layer_attention = attention_weights[layer_idx][sample_idx, head_idx]
    tokens = tokenizer.convert_ids_to_tokens(input_ids[sample_idx])
    
    if len(tokens) > max_tokens:
        tokens = tokens[:max_tokens]
        layer_attention = layer_attention[:max_tokens, :max_tokens]
    
    plt.figure(figsize=figsize)
    
    colors = ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', 
              '#2171b5', '#08519c', '#08306b']
    cmap = LinearSegmentedColormap.from_list('attention', colors, N=100)
    
    sns.heatmap(layer_attention, xticklabels=tokens, yticklabels=tokens,
                cmap=cmap, square=True, cbar_kws={'label': 'Attention Weight'},
                linewidths=0.1)
    
    plt.title(f'Attention Heatmap - Layer {layer_idx}, Head {head_idx}', fontsize=16, pad=20)
    plt.xlabel('Key Tokens', fontsize=12)
    plt.ylabel('Query Tokens', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()


def visualize_multiple_heads(attention_weights, tokenizer, input_ids, layer_idx=-1, 
                            num_heads=16, sample_idx=0, max_tokens=30, cols=4, figsize=(24, 24)):
    """Visualize multiple attention heads in a grid."""
    rows = int(np.ceil(num_heads / cols))
    fig, axes = plt.subplots(rows, cols, figsize=figsize)
    axes = axes.flatten()
    
    tokens = tokenizer.convert_ids_to_tokens(input_ids[sample_idx])
    if len(tokens) > max_tokens:
        tokens = tokens[:max_tokens]
    
    colors = ['#f7fbff', '#deebf7', '#c6dbef', '#9ecae1', '#6baed6', '#4292c6', 
              '#2171b5', '#08519c', '#08306b']
    cmap = LinearSegmentedColormap.from_list('attention', colors, N=100)
    
    for head_idx in range(num_heads):
        if head_idx >= len(axes):
            break
        
        layer_attention = attention_weights[layer_idx][sample_idx, head_idx]
        if len(tokens) < layer_attention.shape[0]:
            layer_attention = layer_attention[:len(tokens), :len(tokens)]
        
        sns.heatmap(layer_attention, xticklabels=tokens, yticklabels=tokens,
                   cmap=cmap, square=True, cbar=False, ax=axes[head_idx],
                   linewidths=0.1)
        
        axes[head_idx].set_title(f'Head {head_idx}', fontsize=10)
        axes[head_idx].tick_params(axis='x', rotation=45, labelsize=7)
        axes[head_idx].tick_params(axis='y', rotation=0, labelsize=7)
    
    plt.suptitle(f'Attention Heads - Layer {layer_idx}', fontsize=18, y=1.02)
    plt.tight_layout()
    plt.show()


def plot_attention_summary(attention_weights, layer_range=None):
    """Plot a summary of attention patterns across layers."""
    num_layers = len(attention_weights)
    if layer_range is None:
        layer_range = range(num_layers)
    
    avg_attentions = []
    for layer_idx in layer_range:
        avg_att = np.mean(attention_weights[layer_idx])
        avg_attentions.append(avg_att)
    
    plt.figure(figsize=(10, 6))
    plt.plot(layer_range, avg_attentions, 'bo-', linewidth=2, markersize=8)
    plt.xlabel('Layer Index', fontsize=12)
    plt.ylabel('Average Attention Weight', fontsize=12)
    plt.title('Average Attention Weight by Layer', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()