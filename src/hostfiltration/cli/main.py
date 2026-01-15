"""
src/Host_Filtration/cli/main.py
Command-line interface for the Host Filtration tool.
"""

import argparse
import sys
import os
from pathlib import Path
from transformers import AutoTokenizer

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from hostfiltration.data.preprocessing import load_data, prepare_datasets
from hostfiltration.models.classifier import (
    setup_device, setup_model, train_model, load_model_from_checkpoint
)
from hostfiltration.models.evaluation import (
    test_model_on_dataset, evaluate_with_trainer,
    visualize_attention_heatmap, visualize_multiple_heads
)


def generate_dataset_command(args):
    """Generate training dataset from genomes."""
    from hostfiltration.data.dataset_generation import main as generate_main
    
    # Override sys.argv to pass arguments to dataset_generation
    sys.argv = [
        'dataset_generation.py',
        '--metadata', args.metadata,
        '--sample-size', str(args.sample_size),
        '--window-size', str(args.window_size),
        '--stride', str(args.stride),
        '--output', args.output,
        '--ratio', args.ratio,
        '--human_genome_path', args.human_genome_path
    ]
    generate_main()


def train_command(args):
    """Train the host filtration model."""
    print("=== Setting up device ===")
    device = setup_device()
    
    print("\n=== Loading data ===")
    train_df, val_df, test_df = load_data(args.data)
    
    print("\n=== Loading tokenizer ===")
    tokenizer = AutoTokenizer.from_pretrained(args.model_name)
    
    print("\n=== Preparing datasets ===")
    train_ds, val_ds, test_ds = prepare_datasets(train_df, val_df, test_df, tokenizer)
    
    print("\n=== Setting up model ===")
    model, config = setup_model(device, model_name=args.model_name)
    
    print("\n=== Training model ===")
    trainer = train_model(
        model=model,
        train_ds=train_ds,
        val_ds=val_ds,
        output_dir=args.output_dir,
        learning_rate=args.learning_rate,
        num_epochs=args.epochs,
        batch_size=args.batch_size
    )
    
    print("\n=== Evaluating on test set ===")
    evaluate_with_trainer(trainer, test_ds)
    
    print(f"\n✓ Training complete! Model saved to {args.output_dir}")


def evaluate_command(args):
    """Evaluate a trained model."""
    print("=== Setting up device ===")
    device = setup_device()
    
    print("\n=== Loading model from checkpoint ===")
    model, config = load_model_from_checkpoint(args.checkpoint, device)
    
    print("\n=== Loading data ===")
    _, _, test_df = load_data(args.data)
    
    print("\n=== Loading tokenizer ===")
    tokenizer = AutoTokenizer.from_pretrained(args.model_name)
    
    print("\n=== Preparing test dataset ===")
    _, _, test_ds = prepare_datasets(test_df, test_df, test_df, tokenizer)
    
    print("\n=== Evaluating model ===")
    acc, f1, cm, attention_weights, failures = test_model_on_dataset(
        model=model,
        ds=test_ds,
        batch_size=args.batch_size,
        threshold=args.threshold
    )
    
    if args.visualize and attention_weights is not None:
        print("\n=== Visualizing attention ===")
        # Get first batch for visualization
        first_batch = test_ds[0]
        input_ids = first_batch['input_ids'].unsqueeze(0)
        
        visualize_attention_heatmap(
            attention_weights=attention_weights,
            tokenizer=tokenizer,
            input_ids=input_ids,
            layer_idx=-1,
            head_idx=0
        )
    
    print(f"\n✓ Evaluation complete!")
    print(f"Accuracy: {acc:.4f}")
    print(f"F1 Score: {f1:.4f}")


def predict_command(args):
    """Predict labels for new sequences."""
    import pandas as pd
    import torch
    import numpy as np
    
    print("=== Setting up device ===")
    device = setup_device()
    
    print("\n=== Loading model from checkpoint ===")
    model, config = load_model_from_checkpoint(args.checkpoint, device)
    model.eval()
    
    print("\n=== Loading tokenizer ===")
    tokenizer = AutoTokenizer.from_pretrained(args.model_name)
    
    # Read input sequences
    if args.input.endswith('.csv'):
        df = pd.read_csv(args.input)
        sequences = df['sequence'].tolist()
    elif args.input.endswith('.fasta') or args.input.endswith('.fa'):
        from Bio import SeqIO
        sequences = [str(record.seq) for record in SeqIO.parse(args.input, 'fasta')]
    else:
        # Single sequence
        sequences = [args.input]
    
    print(f"\n=== Predicting {len(sequences)} sequences ===")
    results = []
    
    with torch.no_grad():
        for seq in sequences:
            # Tokenize
            inputs = tokenizer(seq, padding="max_length", truncation=True, 
                             max_length=512, return_tensors="pt")
            inputs = {k: v.to(device) for k, v in inputs.items()}
            
            # Predict
            outputs = model(**inputs)
            logits = outputs.logits.cpu().numpy()
            probs = np.exp(logits) / np.exp(logits).sum(-1, keepdims=True)
            token_preds = probs.argmax(-1)
            
            # Sequence-level prediction
            human_token_pct = token_preds.sum() / token_preds.shape[1]
            seq_pred = 1 if human_token_pct >= args.threshold else 0
            
            results.append({
                'sequence': seq[:50] + '...' if len(seq) > 50 else seq,
                'prediction': 'Human' if seq_pred == 1 else 'Microbial',
                'human_token_pct': float(human_token_pct),
                'confidence': float(human_token_pct if seq_pred == 1 else 1 - human_token_pct)
            })
    
    # Output results
    results_df = pd.DataFrame(results)
    if args.output:
        results_df.to_csv(args.output, index=False)
        print(f"\n✓ Results saved to {args.output}")
    else:
        print("\n=== Results ===")
        print(results_df.to_string(index=False))


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description=(
            "Host Filtration: ML-based host genome filtering for metagenomics.\n\n"
            "This tool allows you to:\n"
            "  1. Generate a labeled dataset of microbial and human genome subsequences.\n"
            "  2. Train a transformer-based classifier to distinguish human vs microbial sequences.\n"
            "  3. Evaluate a trained model on held-out data.\n"
            "  4. Predict labels for new sequences (CSV, FASTA, or single sequence).\n\n"
            "Example usage:\n"
            "  hostfiltration generate --metadata metadata.tsv --sample-size 10 --window-size 150 --stride 150 --output subsequences_dataset.csv -- ratio 1:1 --human-genome human_genome.fna.qz\n"
            "  hostfiltration train --data dataset.csv --epochs 10\n"
            "  hostfiltration evaluate --checkpoint results/model.ckpt --data test.csv\n"
            "  hostfiltration predict --checkpoint results/model.ckpt --input sequences.fasta"
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Generate dataset command
    generate_parser = subparsers.add_parser(
        'generate',
        help='Generate training dataset from genomes',
        description=(
            "Download microbial genomes listed in a metadata TSV file, extract fixed-length subsequences, "
            "download the human genome, sample random subsequences, and create a labeled dataset CSV.\n\n"
            "Example:\n"
            "  hostfiltration generate --metadata metadata.tsv --sample-size 10 --window-size 150 --stride 150 --output subsequences_dataset.csv -- ratio 1:1"
        )
    )
    generate_parser.add_argument('--metadata', required=True, help='Path to metadata TSV file containing microbial genomes.')
    generate_parser.add_argument('--sample-size', type=int, default=10, help='Number of microbial genomes to sample for dataset creation.')
    generate_parser.add_argument('--window-size', type=int, default=150, help='Length (in bp) of each subsequence.')
    generate_parser.add_argument('--stride', type=int, default=150, help='Step size between subsequences; set < window_size for overlap.')
    generate_parser.add_argument('--output', default='subsequences_dataset.csv', help='Path to save the resulting CSV dataset.')
    generate_parser.add_argument('--ratio', type=str, default='1:1', help='Ratio of human to microbial sequences (e.g., 1:4).')
    generate_parser.add_argument('--human_genome_path', type=str, default=None, help='Path to a local human genome .fna.gz file (optional).')

    # Train command
    train_parser = subparsers.add_parser(
        'train',
        help='Train the host filtration model',
        description=(
            "Train a transformer-based model to classify sequences as human or microbial.\n\n"
            "Example:\n"
            "  hostfiltration train --data dataset.csv --model-name InstaDeepAI/nucleotide-transformer-v2-50m-multi-species"
        )
    )
    train_parser.add_argument('--data', required=True, help='Path to the training dataset CSV.')
    train_parser.add_argument('--model-name', default='InstaDeepAI/nucleotide-transformer-v2-50m-multi-species',
                              help='Name of the pre-trained transformer model to fine-tune.')
    train_parser.add_argument('--output-dir', default='./results', help='Directory to save model checkpoints and outputs.')
    train_parser.add_argument('--learning-rate', type=float, default=1e-5, help='Learning rate for optimizer.')
    train_parser.add_argument('--epochs', type=int, default=5, help='Number of training epochs.')
    train_parser.add_argument('--batch-size', type=int, default=8, help='Batch size for training.')

    # Evaluate command
    eval_parser = subparsers.add_parser(
        'evaluate',
        help='Evaluate a trained model',
        description=(
            "Evaluate a trained model on a test dataset, reporting accuracy, F1 score, and optionally visualize attention weights.\n\n"
            "Example:\n"
            "  hostfiltration evaluate --checkpoint results/model.ckpt --data test.csv --visualize"
        )
    )
    eval_parser.add_argument('--checkpoint', required=True, help='Path to trained model checkpoint.')
    eval_parser.add_argument('--data', required=True, help='Path to test dataset CSV.')
    eval_parser.add_argument('--model-name', default='InstaDeepAI/nucleotide-transformer-v2-50m-multi-species',
                             help='Pre-trained model name used for tokenizer.')
    eval_parser.add_argument('--batch-size', type=int, default=32, help='Batch size for evaluation.')
    eval_parser.add_argument('--threshold', type=float, default=0.1, help='Human token percentage threshold for classification.')
    eval_parser.add_argument('--visualize', action='store_true', help='Visualize attention weights for the first batch.')

    # Predict command
    predict_parser = subparsers.add_parser(
        'predict',
        help='Predict labels for new sequences',
        description=(
            "Predict labels for sequences in CSV, FASTA, or single sequence input using a trained model.\n\n"
            "Example:\n"
            "  hostfiltration predict --checkpoint results/model.ckpt --input sequences.fasta --threshold 0.2"
        )
    )
    predict_parser.add_argument('--checkpoint', required=True, help='Path to trained model checkpoint.')
    predict_parser.add_argument('--input', required=True, help='Input sequences (CSV, FASTA file, or single sequence string).')
    predict_parser.add_argument('--model-name', default='InstaDeepAI/nucleotide-transformer-v2-50m-multi-species',
                                help='Pre-trained model name used for tokenizer.')
    predict_parser.add_argument('--threshold', type=float, default=0.2, help='Human token percentage threshold for classification.')
    predict_parser.add_argument('--output', help='Path to save predictions as CSV. If not provided, prints to console.')

    args = parser.parse_args()

    if args.command == 'generate':
        generate_dataset_command(args)
    elif args.command == 'train':
        train_command(args)
    elif args.command == 'evaluate':
        evaluate_command(args)
    elif args.command == 'predict':
        predict_command(args)
    else:
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()