import os

import configparser
import argparse

import yaml


def parse_ini_file(file):
    config = configparser.ConfigParser(inline_comment_prefixes=[";", "#"])
    with open(file, "r") as f:
        config.read_file(f)

    params = {}
    for sec in config.sections():
        params[sec] = {}
        for key in config[sec].keys():
            val = config[sec][key]
            params[sec][key] = val

    return params


def build_params_spec(values_spec, priors_spec=None, use_cobaya_theory=True, halofit_version=""):
    if priors_spec is None:
        priors_spec = {}

    cobaya_param_spec = {}
    cobaya_cosmology_spec = {}

    for section_name, section_data in values_spec.items():
        for param_name, param_data in section_data.items():
            param_data = list(param_data.split())

            if len(param_data) == 3:
                # Parameter is being sampled

                ref = float(param_data[1])
                # Check if prior is specified
                if section_name in priors_spec and param_name in priors_spec[section_name]:
                    prior_data = list(priors_spec[section_name][param_name].split())
                    prior_distr = prior_data[0].strip()
                    prior_data = prior_data[1:]
                else:
                    prior_distr = "uniform"
                    prior_data = [param_data[0], param_data[2]]
                
                if prior_distr == "gaussian":
                    prior = {
                        "dist": "norm",
                        "loc": float(prior_data[0]),
                        "scale": float(prior_data[1])
                    }
                elif prior_distr == "uniform":
                    prior = {
                        "dist": "uniform",
                        "min": float(prior_data[0]),
                        "max": float(prior_data[1])
                    }
                else:
                    raise RuntimeError(f"Unsupported prior distribution {prior_distr} for {section_name}/{param_name}. Full data: {param_data}")

                cobaya_spec = {
                    "ref": ref,
                    "prior": prior
                }
            elif len(param_data) == 1:
                cobaya_spec = float(param_data[0])
            else:
                raise RuntimeError("Unsupported param spec for {cobaya_name}: {param_data}")

            if use_cobaya_theory:
                # Translate parameter names from CosmoSIS to cobaya conventions
                section_name, param_name = translate_params_to_cobaya_theory(
                    section_name, param_name, halofit_version)

            if section_name == "cosmological_parameters" and use_cobaya_theory:
                # If using the cobaya theory code, move cosmology parameters there
                cobaya_cosmology_spec[param_name] = cobaya_spec
            else:
                cobaya_name = f"{section_name}.{param_name}"
                cobaya_param_spec[cobaya_name] = cobaya_spec
    
    return cobaya_cosmology_spec, cobaya_param_spec


def translate_params_to_cobaya_theory(section_name, param_name, halofit_version=""):
    cobaya_section = section_name
    cobaya_name = param_name

    if section_name == "cosmological_parameters":
        if param_name == "n_s":
            cobaya_name = "ns"
            cobaya_section = "cosmological_parameters"
        elif param_name == "omega_k":
            cobaya_name = "omegak"
            cobaya_section = "cosmological_parameters"
    
    if section_name == "halo_model_parameters":
        # Move HMCode parameters to cobaya theory section for camb
        if param_name == "logT_AGN".lower() and halofit_version == "mead2020_feedback":
            cobaya_name = "HMCode_logT_AGN"
            cobaya_section = "cosmological_parameters"
        elif param_name == "A".lower() and halofit_version in ("mead2015", "mead2016"):
            cobaya_name = "HMCode_A_baryon"
            cobaya_section = "cosmological_parameters"
    
    return cobaya_section, cobaya_name



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pipeline-file", required=True)
    parser.add_argument("--values-file")
    parser.add_argument("--priors-file")
    parser.add_argument("--output-yaml-file", required=True)
    parser.add_argument("--cobaya-output-path", required=True)
    parser.add_argument("--halofit_version", default="mead2020_feedback")
    parser.add_argument("--neutrino_hierarchy", default="normal")
    parser.add_argument("--boltzmann-code", choices=["cobaya", "cosmosis"], default="cobaya")
    parser.add_argument("--sampler", choices=["evaluate", "mcmc"], default="evaluate")
    parser.add_argument("--sampler.mcmc.covmat")

    args = parser.parse_args()

    # Read CosmoSIS main configuration file
    cosmosis_spec = parse_ini_file(args.pipeline_file)

    # Read parameter specification from CosmoSIS values file
    cosmosis_values_file = args.values_file
    if cosmosis_values_file is None:
        print("--values-file not specified, reading from CosmoSIS configuration.")
        cosmosis_values_file = cosmosis_spec["pipeline"]["values"]

    values_spec = parse_ini_file(cosmosis_values_file)

    # Read parameter prior specification from CosmoSIS prior file
    cosmosis_priors_file = args.priors_file
    if cosmosis_priors_file is None:
        print("--priors-file not specified, reading from CosmoSIS configuration.")
        try:
            cosmosis_priors_file = cosmosis_spec["pipeline"]["priors"].strip()
            if cosmosis_priors_file == "":
                cosmosis_priors_file = None
        except KeyError():
            cosmosis_priors_file = None

    priors_spec = {}
    if cosmosis_priors_file is not None:
        priors_spec = parse_ini_file(cosmosis_priors_file)
    else:
        priors_spec = None
    
    use_cobaya_theory = args.boltzmann_code == "cobaya"
    
    camb_settings = {}
    cobaya_theory_spec = {}
    remove_modules = []

    halofit_version = args.halofit_version
    neutrino_hierarchy = args.neutrino_hierarchy

    if use_cobaya_theory:
        remove_modules = ["sample_S8", "sigma8toAs", "one_parameter_hmcode", "camb", "cosmopower", "distances"]

        camb_settings = {}
        if "camb" in cosmosis_spec:
            # Update defaults with values taken from the cosmosis ini
            def cast_options_to_int_float(param_name, value):
                int_options_camb = ["nz", "nz_mid", "nz_background"]
                float_options_camb = ["kmax", "zmax", "zmid", "zmin_background", "zmax_background"]

                if param_name in int_options_camb:
                    return int(value)
                elif param_name in float_options_camb:
                    return float(value)
                else:
                    return value

            camb_options = [
                "kmax",
                
                "zmid",
                "nz_mid",
                "zmax",
                "nz",

                "zmax_background",
                "zmin_background",
                "nz_background",
            ]
            camb_settings = {
                k: cast_options_to_int_float(k, cosmosis_spec["camb"][k])
                for k in camb_options if k in cosmosis_spec["camb"]
            }

            halofit_version = cosmosis_spec["camb"].get("halofit_version", halofit_version)
            neutrino_hierarchy = cosmosis_spec["camb"].get("neutrino_hierarchy", neutrino_hierarchy)
            
        cobaya_theory_spec = {
            "theory": {
                "camb": {
                    "extra_args": {
                        "halofit_version": halofit_version,
                        "neutrino_hierarchy": neutrino_hierarchy,
                    }
                }
            }
        }

    cobaya_cosmology_param_spec, cobaya_param_spec = build_params_spec(
        values_spec, priors_spec, use_cobaya_theory=use_cobaya_theory,
        halofit_version=halofit_version)
    
    if use_cobaya_theory:
        # Translate h0 to H0
        if "h0" in cobaya_cosmology_param_spec:
            cobaya_cosmology_param_spec["H0"] = "lambda h0: h0*100"

        # Transform S8 into sigma8
        if "s_8_input" in cobaya_cosmology_param_spec:
            cobaya_cosmology_param_spec["s_8_input"]["drop"] = True
            cobaya_cosmology_param_spec["sigma8"] = "lambda s_8_input, ombh2, omch2, H0: s_8_input/np.sqrt((ombh2 + omch2)/(H0/100)**2/0.3)"

        cobaya_cosmology_param_spec["omegam"] = {}
        cobaya_cosmology_param_spec["As"] = {}
        cobaya_cosmology_param_spec["s8"] = {"derived": "lambda sigma8, omegam: sigma8*np.sqrt(omegam/0.3)"}

        cobaya_cosmology_param_spec = {"params": cobaya_cosmology_param_spec}
    else:
        cobaya_cosmology_param_spec = {}

    if "extra_output" in cosmosis_spec["pipeline"]:
        derived_param_spec = {}
        for param_name in cosmosis_spec["pipeline"]["extra_output"].split():
            if "#" in param_name:
                # Don't support vector outputs for now
                continue
            section, key = param_name.split("/")
            if use_cobaya_theory and section == "cosmological_parameters":
                # Derived cosmological parameters are computed by cobaya
                continue
            derived_param_spec[f"{section}.{key}"] = {"derived": True}
        
        cobaya_param_spec = {**cobaya_param_spec, **derived_param_spec}

    if args.sampler == "evaluate":
        cobaya_sampler_spec = {"evaluate": None}
    elif args.sampler == "mcmc":
        cobaya_sampler_spec = {"mcmc": None}
        covmat = getattr(args, "sampler.mcmc.covmat")
        if covmat is not None:
            cobaya_sampler_spec["mcmc"] = {"covmat": covmat}

    cobaya_spec = {
        "debug": False,
        "stop_at_error": True,
        "timing": True,
        "output": args.cobaya_output_path,
        "sampler": cobaya_sampler_spec,
        "likelihood": {
            "cosmosis_cobaya_interface.cosmosis_wrapper.CosmoSISWrapperLikelihood": {
                "ini_file": args.pipeline_file,
                "remove_modules": remove_modules,
                "use_cobaya_theory": use_cobaya_theory,
                **camb_settings,
                "params": cobaya_param_spec,
            }
        },
        **cobaya_theory_spec,
        **cobaya_cosmology_param_spec
    }

    # print(yaml.dump(cobaya_spec, sort_keys=False))

    output_dir = os.path.dirname(args.output_yaml_file)
    if output_dir != "" and not os.path.exists(output_dir):
        print(f"Path for output yaml file ({output_dir}) does not exists, creating.")
        os.makedirs(output_dir)
    with open(args.output_yaml_file, "w") as f:
        yaml.dump(cobaya_spec, f, sort_keys=False)
