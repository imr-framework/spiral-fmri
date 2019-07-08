# spiral-fMRI

**Please note:** You need to install Pulseq interpreter on your system before running any Pulseq sequence.

In **PRESTO_pulseq**:
1. Change the system limits if you need to. 
    ```python
    gamma = 42576000  # in Hz/T  %Determined from Pulseq - do not change

    kwargs_for_opts = {'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s', 'grad_dead_time': 10e-6}
    system = Opts(kwargs_for_opts)
    seq = Sequence(system)
    ```
    
    For example, gradient and slew rate limits are not the same for GE and Siemens.
    
    ![GEvsSiemens system limits](/images/rf_g_limits.png)
   
   
2. Run the code . It will give you the .seq file that will be played on the scanner.
    ```Python
    seq.write("spiral_PRESTO.seq")
    ```
3. Go to the MR room and play it!

