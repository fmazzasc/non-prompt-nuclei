if [ ! -d "results" ]; then
  mkdir results
fi

# python3 prepare_deut_data.py --config-file configs/config_d_thnsparse.yaml
# python3 dca_analysis.py --config-file configs/config_d_thnsparse.yaml
# python3 prepare_deut_data.py --config-file configs/config_d_tree.yaml
# python3 dca_analysis.py --config-file configs/config_d_tree.yaml

python3 prepare_he3_data.py --config-file configs/config_he3.yaml
python3 dca_analysis.py --config-file configs/config_he3.yaml