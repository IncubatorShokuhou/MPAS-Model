#!/usr/bin/env bash
# Jablonowski-Williamson baroclinic wave (init case 2) baseline.
#
# Two meshes:
#   uniform   x1.40962   120 km, 40962 cells
#             grid from the NCAR idealized-case tarball (v7.0), which is the
#             same quasi-uniform mesh NCAR publishes as the 120 km download
#   variable  x5.30210   240 km -> 48 km (factor 5), 30210 cells
#             hosted with the other NCAR meshes but not listed on the mesh table;
#             this is the mesh the MPAS tutorial calls x5.30210
#
# Dynamics knobs follow the v7.0 jw_baroclinic_wave namelist, with one v8 change:
# config_epssm was removed from the integration (it only warns). Constant
# off-centering of 0.1 is reproduced with config_epssm_minimum = config_epssm_maximum = 0.1.
#
# Usage:
#   ./docs/hts/scripts/run_jw_baseline.sh [uniform|variable|both] [duration]
# Example:
#   ./docs/hts/scripts/run_jw_baseline.sh both 10_00:00:00
#
# Environment:
#   MPAS_SRC    model root containing atmosphere_model and init_atmosphere_model
#   MPAS_DATA   directory with meshes/ and jw/jw_baroclinic_wave/
#   MPAS_RUNS   where case directories are created
#   MPAS_NP     MPI ranks (default 4). A graph.info.part.${MPAS_NP} must exist.
#   MPAS_OUTPUT_INTERVAL  history interval (default 1_00:00:00, one file per day)

set -euo pipefail

WHICH=${1:-both}
DURATION=${2:-1_00:00:00}
OUTPUT_INTERVAL=${MPAS_OUTPUT_INTERVAL:-1_00:00:00}
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MPAS_SRC=${MPAS_SRC:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
MPAS_DATA=${MPAS_DATA:-"$HOME/mpas-data"}
MPAS_RUNS=${MPAS_RUNS:-"$HOME/mpas-runs"}
MPAS_NP=${MPAS_NP:-4}

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
# Shared-memory BTL can fail inside some VMs.
export OMPI_MCA_btl_vader_single_copy_mechanism=${OMPI_MCA_btl_vader_single_copy_mechanism:-none}

if [ ! -x "$MPAS_SRC/atmosphere_model" ] || [ ! -x "$MPAS_SRC/init_atmosphere_model" ]; then
    echo "ERROR: build atmosphere_model and init_atmosphere_model in $MPAS_SRC first." >&2
    exit 1
fi

write_init_namelist() {
    local dest=$1 prefix=$2
    cat > "$dest/namelist.init_atmosphere" <<EOF
&nhyd_model
   config_start_time = '0000-01-01_00:00:00'
   config_init_case = 2
/
&dimensions
   config_nvertlevels = 26
/
&decomposition
   config_block_decomp_file_prefix = '${prefix}.graph.info.part.'
/
EOF
}

write_init_streams() {
    local dest=$1
    cat > "$dest/streams.init_atmosphere" <<'EOF'
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="grid.nc"
                  input_interval="initial_only" />

<immutable_stream name="output"
                  type="output"
                  filename_template="init.nc"
                  packages="initial_conds"
                  output_interval="initial_only" />
</streams>
EOF
}

write_atm_namelist() {
    local dest=$1 prefix=$2 dt=$3 len_disp=$4
    cat > "$dest/namelist.atmosphere" <<EOF
&nhyd_model
    config_dt = ${dt}
    config_start_time = '0000-01-01_00:00:00'
    config_run_duration = '${DURATION}'
    config_split_dynamics_transport = false
    config_number_of_sub_steps = 6
    config_dynamics_split_steps = 1
    config_h_mom_eddy_visc2    = 0.0
    config_h_mom_eddy_visc4    = 0.0
    config_v_mom_eddy_visc2    = 0.0
    config_h_theta_eddy_visc2  = 0.0
    config_h_theta_eddy_visc4  = 0.0
    config_v_theta_eddy_visc2  = 0.0
    config_horiz_mixing        = '2d_smagorinsky'
    config_len_disp            = ${len_disp}
    config_u_vadv_order      = 3
    config_w_vadv_order      = 3
    config_theta_vadv_order  = 3
    config_scalar_vadv_order = 3
    config_theta_adv_order   = 3
    config_scalar_adv_order  = 3
    config_scalar_advection  = false
    config_positive_definite = false
    config_coef_3rd_order = 1.0
    config_monotonic = false
    config_smdiv = 0.1
/

&damping
    config_zd = 22000.0
    config_xnutr = 0.0
    config_epssm_minimum = 0.1
    config_epssm_maximum = 0.1
    config_epssm_transition_bottom_z = 30000.
    config_epssm_transition_top_z = 50000.
/

&limited_area
    config_apply_lbcs = false
/

&io
    config_pio_num_iotasks = 0
    config_pio_stride = 1
/

&decomposition
    config_block_decomp_file_prefix = '${prefix}.graph.info.part.'
/

&restart
    config_do_restart = false
/

&printout
    config_print_global_minmax_vel = true
    config_print_global_minmax_sca = false
/

&physics
    config_sst_update = false
    config_sstdiurn_update = false
    config_deepsoiltemp_update = false
    config_radtlw_interval = '00:30:00'
    config_radtsw_interval = '00:30:00'
    config_bucket_update = 'none'
    config_physics_suite = 'none'
/
EOF
}

write_atm_streams() {
    local dest=$1
    # One history file. output_interval equals the run length when DURATION
    # uses the DDD_hh:mm:ss form; also accept hh:mm:ss.
    cat > "$dest/streams.atmosphere" <<EOF
<streams>
<immutable_stream name="input"
                  type="input"
                  filename_template="init.nc"
                  input_interval="initial_only"/>

<immutable_stream name="restart"
                  type="input;output"
                  filename_template="restart.nc"
                  input_interval="initial_only"
                  output_interval="none"/>

<stream name="output"
        type="output"
        filename_template="history.nc"
        filename_interval="none"
        output_interval="${OUTPUT_INTERVAL}">
    <var name="u"/>
    <var name="w"/>
    <var name="rho"/>
    <var name="theta"/>
    <var name="pressure"/>
    <var name="surface_pressure"/>
    <var name="divergence"/>
    <var name="vorticity"/>
    <var name="areaCell"/>
    <var name="zgrid"/>
    <var name="meshDensity"/>
    <var name="dcEdge"/>
    <var name="latCell"/>
    <var name="lonCell"/>
</stream>

<immutable_stream name="iau"
                  type="input"
                  filename_template="iau_not_used.nc"
                  packages="iau"
                  input_interval="none" />

<immutable_stream name="lbc_in"
                  type="input"
                  filename_template="lbc_not_used.nc"
                  packages="limited_area"
                  input_interval="none" />
</streams>
EOF
}

prepare_case() {
    local name=$1 grid=$2 graphdir=$3 prefix=$4 dt=$5 len_disp=$6
    local dest="$MPAS_RUNS/$name"
    rm -rf "$dest"
    mkdir -p "$dest"
    ln -sfn "$MPAS_SRC/init_atmosphere_model" "$dest/init_atmosphere_model"
    ln -sfn "$MPAS_SRC/atmosphere_model" "$dest/atmosphere_model"
    ln -sfn "$grid" "$dest/grid.nc"
    local part
    for part in "$graphdir/${prefix}.graph.info.part."*; do
        ln -sfn "$part" "$dest/$(basename "$part")"
    done
    if [ ! -e "$dest/${prefix}.graph.info.part.${MPAS_NP}" ]; then
        echo "ERROR: missing ${prefix}.graph.info.part.${MPAS_NP} in $graphdir" >&2
        exit 1
    fi
    write_init_namelist "$dest" "$prefix"
    write_init_streams "$dest"
    write_atm_namelist "$dest" "$prefix" "$dt" "$len_disp"
    write_atm_streams "$dest"
    echo "$dest"
}

run_case() {
    local name=$1
    local dest="$MPAS_RUNS/$name"
    echo "========== INIT $name =========="
    (
        cd "$dest"
        /usr/bin/time -f "INIT_WALL_SEC %e\nINIT_MAXRSS_KB %M" -o time.init.txt \
            mpirun -np 1 ./init_atmosphere_model
    )
    if ! grep -q "Finished running the init_atmosphere core" "$dest/log.init_atmosphere.0000.out"; then
        echo "ERROR: init failed for $name. See $dest/log.init_atmosphere.0000.out" >&2
        tail -40 "$dest/log.init_atmosphere.0000.out" || true
        exit 1
    fi
    echo "========== MODEL $name  np=$MPAS_NP  duration=$DURATION =========="
    (
        cd "$dest"
        /usr/bin/time -f "MODEL_WALL_SEC %e\nMODEL_USER_SEC %U\nMODEL_SYS_SEC %S\nMODEL_MAXRSS_KB %M" \
            -o time.model.txt \
            mpirun -np "$MPAS_NP" ./atmosphere_model
    )
    if ! grep -q "Finished running the atmosphere core" "$dest/log.atmosphere.0000.out"; then
        echo "ERROR: atmosphere_model failed for $name" >&2
        tail -50 "$dest/log.atmosphere.0000.out" || true
        exit 1
    fi
    echo "----- $name timings -----"
    cat "$dest/time.init.txt" "$dest/time.model.txt"
}

case "$WHICH" in
    uniform|both)
        prepare_case jw_uniform_120km \
            "$MPAS_DATA/jw/jw_baroclinic_wave/x1.40962.grid.nc" \
            "$MPAS_DATA/jw/jw_baroclinic_wave" \
            x1.40962 \
            450.0 \
            120000.0
        ;;
esac
case "$WHICH" in
    variable|both)
        # dt and len_disp follow the finest nominal spacing (48 km).
        # 120 km archive uses dt=450 s and len_disp=120000 m.
        # 450 * (48/120) = 180 s; len_disp = nominalMinDc * earth radius = 48000 m.
        prepare_case jw_variable_240_48km \
            "$MPAS_DATA/meshes/x5.30210.grid.nc" \
            "$MPAS_DATA/meshes" \
            x5.30210 \
            180.0 \
            48000.0
        ;;
esac

case "$WHICH" in
    uniform) run_case jw_uniform_120km ;;
    variable) run_case jw_variable_240_48km ;;
    both)
        run_case jw_uniform_120km
        run_case jw_variable_240_48km
        ;;
    *)
        echo "Usage: $0 [uniform|variable|both] [duration]" >&2
        exit 1
        ;;
esac
