# spiral-fMRI
This project presents a spiral-based functional Magnetic Resonance Imaging (fMRI) sequence. It allows dynamic T2* weighed contrast with a high temporal resolution and shorter total acquisition time.5 

This Pulseq version is a "translation" of the original one based in TOPPE and it demonstrates the feasibility of porting sequences between these two open-source and vendor-independent frameworks. It has also been implemented in a cross-vendor manner.

Check the [Wiki](https://github.com/imr-framework/spiral-fmri/wiki) for more information.
## Dependencies
## Demo
## Contributing and Community guidelines
## References

**Please note:** You need to install Pulseq interpreter on your system before running any Pulseq sequence. Visit their [GitHub page](http://pulseq.github.io/) for more information

In **PRESTO_pulseq.py**:
1. Change the system limits if you need to. We have tested our implementation with these limits:
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

Check the [Wiki](https://github.com/imr-framework/spiral-fmri/wiki) for more information.
