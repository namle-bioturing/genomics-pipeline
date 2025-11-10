```bash
/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wgs_samples_test.tsv wgs /mnt/disk1/namle/run/giab/output/wgs_pcr_free

/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wgs_samples.tsv wgs /mnt/disk1/namle/run/giab/output/wgs_pcr_free
```

```bash
# Start a tmux session
tmux new-session -s giab-wes-idt

# Attach later
tmux attach -t giab-wes-idt

# Remove 
tmux kill-session -t giab-wes-idt

# View logs
tmux capture-pane -pt giab-wes-idt -S - | less
```

```bash
/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wes_samples_idt_50x.tsv wes /mnt/disk1/namle/run/giab/output/wes_idt_50x

/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wes_samples_idt_75x.tsv wes /mnt/disk1/namle/run/giab/output/wes_idt_75x

/bin/bash /mnt/disk1/namle/run/giab/run_samples.sh /mnt/disk1/namle/run/giab/wes_samples_idt_100x.tsv wes /mnt/disk1/namle/run/giab/output/wes_idt_100x
```