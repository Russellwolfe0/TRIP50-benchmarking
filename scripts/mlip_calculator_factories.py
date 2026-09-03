"""Lazy-loaded ASE calculator providers for the TRIP50 MLIP driver.

Every backend is imported only when selected, allowing one driver to live in
all of the otherwise separate model environments.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable


@dataclass
class CalculatorProvider:
    """Supply either one shared calculator or one calculator per structure."""

    calculator: Any = None
    per_structure: Callable[[Any], Any] | None = None

    def for_atoms(self, atoms: Any) -> Any:
        return self.per_structure(atoms) if self.per_structure else self.calculator


def fairchem(model: str, device: str, **kwargs: Any) -> CalculatorProvider:
    try:
        from fairchem.core import FAIRChemCalculator
        try:
            from fairchem.core.calculate import pretrained_mlip
        except ImportError:
            from fairchem.core import pretrained_mlip
    except ImportError as error:
        raise RuntimeError("Install/activate fairchem-core for the FairChem backend") from error
    task_name = kwargs.pop("task_name", "omol")
    predictor = pretrained_mlip.get_predict_unit(model, device=device, **kwargs)
    return CalculatorProvider(FAIRChemCalculator(predictor, task_name=task_name))


def mace(model: str, device: str, **kwargs: Any) -> CalculatorProvider:
    try:
        from mace import calculators
    except ImportError as error:
        raise RuntimeError("Install/activate mace-torch for the MACE backend") from error
    # MACE exposes OMOL and Polar through separate factory functions. Infer the
    # Polar family from its unambiguous checkpoint names, while retaining an
    # explicit ``family`` override for local/custom checkpoints.
    family = kwargs.pop("family", None)
    if family is None:
        family = "polar" if model.lower().startswith("polar-") else "omol"
    try:
        factory = getattr(calculators, f"mace_{family}")
    except AttributeError as error:
        raise ValueError(f"Unknown MACE family: {family!r}") from error
    kwargs.setdefault("default_dtype", "float64")
    return CalculatorProvider(factory(model=model, device=device, **kwargs))


def aimnet(model: str, device: str, **kwargs: Any) -> CalculatorProvider:
    try:
        from aimnet.calculators import AIMNet2ASE, AIMNet2Calculator
    except ImportError as error:
        raise RuntimeError("Install/activate AIMNet2 for the AIMNet backend") from error
    # Current AIMNet loaders select their own torch device. Accepting device in
    # the common CLI still keeps invocations uniform across backends.
    del device
    base = AIMNet2Calculator(model, **kwargs)

    def for_atoms(atoms: Any) -> Any:
        return AIMNet2ASE(
            base,
            charge=atoms.info.get("charge", 0),
            mult=int(atoms.info.get("spin", 1)),
        )

    return CalculatorProvider(per_structure=for_atoms)


def ani(model: str, device: str, **kwargs: Any) -> CalculatorProvider:
    try:
        import torchani
    except ImportError as error:
        raise RuntimeError("Install/activate TorchANI for the ANI backend") from error
    try:
        model_factory = getattr(torchani.models, model)
    except AttributeError as error:
        raise ValueError(f"Unknown TorchANI model: {model!r} (for example ANI2x)") from error
    return CalculatorProvider(model_factory(**kwargs).to(device).ase())


def orb(model: str, device: str, **kwargs: Any) -> CalculatorProvider:
    try:
        from orb_models.forcefield import pretrained
        from orb_models.forcefield.inference.calculator import ORBCalculator
    except ImportError as error:
        raise RuntimeError("Install/activate orb-models for the ORB backend") from error
    precision = kwargs.pop("precision", "float32-high")
    try:
        model_factory = getattr(pretrained, model)
    except AttributeError as error:
        raise ValueError(f"Unknown ORB pretrained model: {model!r}") from error
    forcefield, atoms_adapter = model_factory(device=device, precision=precision, **kwargs)
    return CalculatorProvider(
        ORBCalculator(forcefield, atoms_adapter=atoms_adapter, device=device)
    )


def mattersim(model: str, device: str, **kwargs: Any) -> CalculatorProvider:
    try:
        from mattersim.forcefield import MatterSimCalculator
    except ImportError as error:
        raise RuntimeError("Install/activate MatterSim for the MatterSim backend") from error
    # The 1M checkpoint is MatterSim's default. For 5M, --checkpoint should be
    # used because installations may store it in different locations.
    checkpoint = kwargs.pop("checkpoint", None)
    if model.lower() not in {"1m", "mattersim-1m", "default"} and checkpoint is None:
        checkpoint = model
    if checkpoint is not None:
        kwargs["load_path"] = checkpoint
    return CalculatorProvider(MatterSimCalculator(device=device, **kwargs))


BACKENDS = {
    "aimnet": aimnet,
    "ani": ani,
    "fairchem": fairchem,
    "mace": mace,
    "mattersim": mattersim,
    "orb": orb,
}


def build_provider(backend: str, model: str, device: str, **kwargs: Any) -> CalculatorProvider:
    """Construct a provider from the stable backend/model interface."""
    try:
        factory = BACKENDS[backend.lower()]
    except KeyError as error:
        choices = ", ".join(sorted(BACKENDS))
        raise ValueError(f"Unknown backend {backend!r}; choose one of: {choices}") from error
    return factory(model=model, device=device, **kwargs)
