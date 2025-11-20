"""
Model training and evaluation module for Host Filtration.
"""

# Model training
from .classifier import (
    setup_device,
    setup_model,
    load_model_from_checkpoint,
    LossLoggingCallback,
    compute_metrics,
    train_model
)

# Model evaluation
from .evaluation import (
    get_sequence_labels,
    test_model_on_dataset,
    evaluate_with_trainer,
    visualize_attention_heatmap,
    visualize_multiple_heads,
    plot_attention_summary
)

__all__ = [
    # model.py
    "setup_device",
    "setup_model",
    "load_model_from_checkpoint",
    "LossLoggingCallback",
    "compute_metrics",
    "train_model",
    # evaluation.py
    "get_sequence_labels",
    "test_model_on_dataset",
    "evaluate_with_trainer",
    "visualize_attention_heatmap",
    "visualize_multiple_heads",
    "plot_attention_summary"
]
