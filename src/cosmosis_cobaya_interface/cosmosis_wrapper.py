import numpy as np

from cobaya import Likelihood

from cosmosis.runtime.config import Inifile
from cosmosis.runtime.pipeline import Pipeline
from cosmosis.datablock.cosmosis_py import block, section_names


class CosmoSISWrapperLikelihood(Likelihood):
    ini_file: str
    remove_modules: list[str] = []
    use_cobaya_theory: bool = True
    dump_datablock_path: str | None = None

    kmin_extrapolate: float = 1e-5
    kmax: float = 10.0
    kmax_extrapolate: float | None = None
    nk: int = 200

    zmin: float = 0.0
    zmax: float = 3.01
    zmid: float | None = None

    nz: int = 150
    nz_mid: int = 100

    zmin_background: float | None = None
    zmax_background: float | None = None
    nz_background: float | None = None

    def initialize(self):
        ini = Inifile(filename=self.ini_file)
        ini.set("pipeline", "quiet", "T")
        ini.set("pipeline", "timing", "F")

        modules = ini.get("pipeline", "modules")
        modules = modules.strip().split(" ")
        for module in self.remove_modules:
            if module in modules:
                modules.remove(module)
                print(f"Removing module {module}")

        ini.set("pipeline", "modules", " ".join(modules))

        if self.use_cobaya_theory:
            if self.zmid is not None:
                self.z_pk = np.concatenate(
                    (
                        np.linspace(self.zmin, self.zmid, self.nz_mid, endpoint=False),
                        np.linspace(self.zmid, self.zmax, self.nz - self.nz_mid),
                    )
                )[::-1]
            else:
                self.z_pk = np.linspace(self.zmin, self.zmax, self.nz)[::-1]

            self.zmin_background = (
                self.zmin if self.zmin_background is None else self.zmin_background
            )
            self.zmax_background = (
                self.zmax if self.zmax_background is None else self.zmax_background
            )
            self.nz_background = (
                self.nz if self.nz_background is None else self.nz_background
            )
            self.z_background = np.linspace(
                self.zmin_background, self.zmax_background, self.nz_background
            )

            self.kmax_extrapolate = (
                self.kmax if self.kmax_extrapolate is None else self.kmax_extrapolate
            )

        self.pipeline = Pipeline(ini)

        self.pipeline.setup()

    def get_requirements(self):
        if self.use_cobaya_theory:
            return {
                "H0": None,
                "omegam": None,
                "omch2": None,
                "ombh2": None,
                "Pk_interpolator": {
                    "z": self.z_pk,
                    "vars_pairs": [["delta_tot", "delta_tot"]],
                    "k_max": self.kmax,
                    "nonlinear": (False, True),
                },
                "comoving_radial_distance": {"z": self.z_background},
            }
        else:
            return {}

    def logp(self, **params):
        data = block.DataBlock()
        # Get cosmosis params
        for name, value in params.items():
            if "." in name:
                section, key = name.split(".")
                data[section, key] = value
            # elif "_derived" not in name:
            #     data[section_names.cosmological_parameters, name] = value

        if self.use_cobaya_theory:
            h = params["H0"] / 100

            data[section_names.cosmological_parameters, "h"] = h
            data[section_names.cosmological_parameters, "h0"] = h
            data[section_names.cosmological_parameters, "omega_m"] = params["omegam"]

            # Get distances
            data[section_names.distances, "z"] = self.z_background
            data[section_names.distances, "a"] = 1 / (1 + self.z_background)
            data[
                section_names.distances, "d_m"
            ] = self.provider.get_comoving_radial_distance(self.z_background)

            # Get matter power spectra
            kmax = max(self.kmax, self.kmax_extrapolate)
            k = np.geomspace(self.kmin_extrapolate, kmax, self.nk) * h
            z = self.z_pk[::-1]
            pk_lin = self.provider.get_Pk_interpolator(
                extrap_kmin=self.kmin_extrapolate * h,
                extrap_kmax=self.kmax_extrapolate * h,
                nonlinear=False,
            ).P(z, k)
            pk_nonlin = self.provider.get_Pk_interpolator(
                extrap_kmin=self.kmin_extrapolate * h,
                extrap_kmax=self.kmax_extrapolate * h,
                nonlinear=True,
            ).P(z, k)

            data.put_grid(
                section=section_names.matter_power_lin,
                name_x="z",
                x=z,
                name_y="k_h",
                y=k / h,
                name_z="p_k",
                z=pk_lin * h**3,
            )
            data.put_grid(
                section=section_names.matter_power_nl,
                name_x="z",
                x=z,
                name_y="k_h",
                y=k / h,
                name_z="p_k",
                z=pk_nonlin * h**3,
            )

        # data[section_names.growth_parameters, ""] = None

        self.pipeline.run(data)

        if self.dump_datablock_path is not None:
            data.save_to_directory(self.dump_datablock_path)

        # for s in data.sections():
        #     print(data.keys(s))
        # print(data[section_names.likelihoods, "log_like"])
        return data[section_names.likelihoods, "loglike_like"]
