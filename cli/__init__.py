# Add metadata or import CLI components

"""
CLI entry points for SCOOTI tools:
- flux_predict_matlab
- infer_objective
- flux_sampler
"""
from .flux_predict_matlab import app as flux_predict_app
from .infer_objective import app as infer_objective_app
from .flux_sampler import app as flux_sampler_app

