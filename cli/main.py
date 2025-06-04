# cli/main.py
import typer
from . import flux_predict_matlab, infer_objective, flux_sampler

app = typer.Typer()
app.add_typer(flux_predict_matlab.app, name="predict")
app.add_typer(infer_objective.app, name="infer")
app.add_typer(flux_sampler.app, name="sample")

if __name__ == "__main__":
    app()

