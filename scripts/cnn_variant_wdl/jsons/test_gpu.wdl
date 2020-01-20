version 1.0

workflow TestGPU{
    File monitoring_script="gs://fc-secure-cbdbcc1d-6e0c-421b-9629-cb8b7df992cd/cromwell_monitoring_script.sh"
    call RunGPUTest{input:
    monitoring_script=monitoring_script}
}

task RunGPUTest{
    input{
    File monitoring_script
    }
    command {
        bash ~{monitoring_script} > monitoring.log &
        sleep 150

        nvidia-smi -q -d UTILIZATION -i 0 | \
                    grep -m2 -e Memory -e Gpu | \
                    sed 's/Memory    /GPU-Memory/';
    }
    runtime {
      docker: "samfriedman/gpu:latest"
      gpuType: "nvidia-tesla-k80" # This will require PAPI v2 and CUDA on VM
      gpuCount: 1
      zones: ["us-central1-c"]
      memory: 1000 + " MB"
      disks: "local-disk 400 SSD"
      bootDiskSizeGb: "16"
    }
    output{
        File monitoring_log="monitoring.log"
    }
}