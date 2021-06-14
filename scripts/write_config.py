import os.path
import yaml

import plasmid_rep


config = {
    'ebv': {
        's_phase_duplication_prob': 0.88,
        'plasmid_repulsion_attraction': 0.24,
        'cluster_jostling': 1.0,
        'cluster_crp_alpha': 101,
        'cluster_jostling_s_vs_g1': 1.0,
        'average_cell_replication_prob': 1.0,
        'signal_selective_disadvantage_on_cell_replication': 0,
        'signal_selective_disadvantage_on_cell_replication_squared': 0,
        'positive_selection_coefficient': 0.1,
        'negative_cluster_selection_coefficient': 0.0
    },

    'kshv': {
        's_phase_duplication_prob': 1.00,  # Equivalent to 0.91 with negative cluster selection coefficient at 0.08
        'plasmid_repulsion_attraction': 1.0,
        'cluster_jostling': 0.8,
        'cluster_crp_alpha': 0.5,
        'cluster_jostling_s_vs_g1': 1.0,
        'average_cell_replication_prob': 1.0,
        'signal_selective_disadvantage_on_cell_replication': 0,
        'signal_selective_disadvantage_on_cell_replication_squared': 0,
        'positive_selection_coefficient': 0.1,
        'negative_cluster_selection_coefficient': 0.07
    },

    'base_path': os.path.join(os.path.expanduser('~'), 'Desktop')
}


def write_config(cfg: dict):
    lib_base, _ = os.path.split(os.path.abspath(plasmid_rep.__file__))
    cfg_path = os.path.join(lib_base, 'config.yaml')
    with open(cfg_path, 'w') as fo:
        yaml.dump(cfg, fo)


if __name__ == '__main__':
    write_config(config)
