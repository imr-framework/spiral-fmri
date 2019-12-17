<h1 align="center"> Spiral fMRI </h1> <br>

This project presents a spiral-based functional Magnetic Resonance Imaging (fMRI) sequence. It allows dynamic T2* weighed contrast with a high temporal resolution and shorter total acquisition time [[1]](#References).

This Pulseq [[2]](#References)[[3]](#References)  version is a "translation" of the original one based in TOPPE [[4]](#References) and it demonstrates the feasibility of porting sequences between these two open-source and vendor-independent frameworks. It has also been implemented in a cross-vendor manner.

Check the [Wiki](https://github.com/imr-framework/spiral-fmri/wiki) for more information.
## Dependencies
1. PyPulseq>=1.2.1
1. numpy>=1.16.3
2. matplotlib>=3.0.3
## Quick Demo
**Please note:** You need to install Pulseq interpreter on your system before running any Pulseq sequence. Visit their [GitHub page](http://pulseq.github.io/) and [specifications](http://pulseq.github.io/specification.pdf) document for more information.
1. Create the environment: 

    After cloning the respository, create a virtual environment, activate it and install the specified [dependencies](#Dependencies).
     ```code
    # Windows
    py -m pip install --user virtualenv # install vitual environment package
    py -m venv <virtual_environment_name> # create a virtual environment in the project directory
    .\<virtual_environment_name>\Scripts\activate # activate it
    pip install -r requirements.txt # install the required packages
    
    # Mac/Linux
    python3 -m pip install --user virtualenv 
    python3 -m venv <virtual_environment_name> 
    source <virtual_environment_name>/bin/activate 
    pip install -r requirements.txt 
    ```
    
    Go to [PyPulseq's GitHub page](https://github.com/imr-framework/pypulseq) to learn more about this tool for pulse sequence design.
2. Open **getparams.py** and check the acquisition parameters

    Note that the sequence implementation has been valited with the default parameters. A different combination of values would affect the resulting images. 
    
    Change the system limits if you need to. We have tested our implementation with these limits:
    ```python
    # system limits
    gamma = 42576000 # Hz/T
    sysSiemens = {'gamma': gamma, 'max_grad': 32, 'grad_unit': 'mT/m', 'max_slew': 130, 'slew_unit': 'T/m/s',
                 'grad_raster_time': 10e-6, 'rf_raster_time': 1e-6}
    ```
    
    For example, gradient and slew rate limits are not the same for GE and Siemens.
    
    ![GEvsSiemens system limits](/images/rf_g_limits.png)
    
3. Run **spiral_fmri_pypulseq.py**. It will give you the .seq file that will be played on the scanner.
    ```Python
    seq.write("spiral_PRESTO.seq")
    ```
4. Go to the MR room and play it!
    
## Contributing and Community guidelines
This repository adheres to a code of conduct adapted from the [Contributor Covenant](https://www.contributor-covenant.org/) code of conduct. Contributing guidelines can be found [here](https://github.com/imr-framework/spiral-fmri/blob/master/CONTRIBUTING.md)
## References
1. van Gelderen, P., Duyn, J. H., Ramsey, N. F., Liu, G., & Moonen, C. T. (2012). The PRESTO technique for fMRI. NeuroImage, 62(2), 676–681. doi:10.1016/j.neuroimage.2012.01.017

2.	Ravi, K. S., Potdar, S., Poojar, P., Reddy, A. K., Kroboth, S., Nielsen, J. F., Zaitsev, M., Venkatesan, R., Geethanath, S. (2018). Pulseq-Graphical Programming Interface: Open source visual environment for prototyping pulse sequences and integrated magnetic resonance imaging algorithm development. Magnetic Resonance Imaging, 52, 9-15. doi:10.1016/j.mri.2018.03.008.

3. Layton, K. J., Kroboth, S., Jia, F., Littin, S., Yu, H., Leupold, J., Nielsen, J., Stöcker, T. and Zaitsev, M. (2017), Pulseq: A rapid and hardware‐independent pulse sequence prototyping framework. Magn. Reson. Med., 77: 1544-1552. doi:10.1002/mrm.26235

4. Nielsen, J. F., & Noll, D. C. (2018). TOPPE: A framework for rapid prototyping of MR pulse sequences. Magnetic resonance in medicine, 79(6), 3128–3134. doi:10.1002/mrm.26990
