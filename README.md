# spiral-fMRI
This project presents a spiral-based functional Magnetic Resonance Imaging (fMRI) sequence. It allows dynamic T2* weighed contrast with a high temporal resolution and shorter total acquisition time. [1](##References) 

This Pulseq[2][3] version is a "translation" of the original one based in TOPPE[4] and it demonstrates the feasibility of porting sequences between these two open-source and vendor-independent frameworks. It has also been implemented in a cross-vendor manner.

Check the [Wiki](https://github.com/imr-framework/spiral-fmri/wiki) for more information.
## Dependencies
## Demo
## Contributing and Community guidelines
## References
1. van Gelderen, P., Duyn, J. H., Ramsey, N. F., Liu, G., & Moonen, C. T. (2012). The PRESTO technique for fMRI. NeuroImage, 62(2), 676–681. doi:10.1016/j.neuroimage.2012.01.017

2.	Ravi, K. S., Potdar, S., Poojar, P., Reddy, A. K., Kroboth, S., Nielsen, J. F., Zaitsev, M., Venkatesan, R., Geethanath, S. (2018). Pulseq-Graphical Programming Interface: Open source visual environment for prototyping pulse sequences and integrated magnetic resonance imaging algorithm development. Magnetic Resonance Imaging, 52, 9-15. doi:10.1016/j.mri.2018.03.008.

3. Layton, K. J., Kroboth, S., Jia, F., Littin, S., Yu, H., Leupold, J., Nielsen, J., Stöcker, T. and Zaitsev, M. (2017), Pulseq: A rapid and hardware‐independent pulse sequence prototyping framework. Magn. Reson. Med., 77: 1544-1552. doi:10.1002/mrm.26235

4. Nielsen, J. F., & Noll, D. C. (2018). TOPPE: A framework for rapid prototyping of MR pulse sequences. Magnetic resonance in medicine, 79(6), 3128–3134. doi:10.1002/mrm.26990

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
