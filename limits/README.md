# pyhf-tutorial

Here we have three notebooks that are independent of each other, which will teach you different things about statistical inference with [pyhf](https://pyhf.readthedocs.io/en/v0.7.2/#) and convenience packages such as [cabinetry](https://cabinetry.readthedocs.io/en/latest/index.html).

The different notebooks include the following:
* [**01-histogram-fits**](./01-histogram-fits.ipynb): This introduces the basics of [pyhf](https://pyhf.readthedocs.io/en/v0.7.2/#), how the model building works and how to include uncertainties. **We encourage you to start with this notebook.**
* [**02-B2Kpi**](./02-B2Kpi.ipynb): This notebook is a realistic example of how to build a statistical model for the $B^+ \to K^+ \pi^0$ decay. We use reconstructed MC in combination with [cabinetry](https://cabinetry.readthedocs.io/en/latest/index.html) to build our [pyhf](https://pyhf.readthedocs.io/en/v0.7.2/#) model. Additionally, tracking efficiency and PID systematics are included in the model. 
* [**03-hypothesis-testing**](./03-hypothesis-testing.ipynb): This notebook goes beyond the basics in [pyhf](https://pyhf.readthedocs.io/en/v0.7.2/#), and introduces advanced methods of statistical inference, such as hypothesis testing on a very simple model.

## References
### pyhf
* [Documentation](https://pyhf.readthedocs.io/en/v0.7.2/#)
* [Overview slides](https://indico.belle2.org/event/8470/contributions/55827/attachments/21257/31463/pyhf.pdf)
* [`pyhf` tutorial](https://pyhf.github.io/pyhf-tutorial/introduction.html)
### HistFactory and asymptotic formulae
* [HistFactory paper](https://cds.cern.ch/record/1456844/files/CERN-OPEN-2012-016.pdf)
* [Asymptotic formulae for likelihood-based tests of new physics](https://arxiv.org/pdf/1007.1727.pdf)

### Cabinetry
* [Documentation](https://cabinetry.readthedocs.io/en/latest/index.html)
* [`Cabinetry` tutorial](https://github.com/cabinetry/cabinetry-tutorials/blob/master/example.ipynb)
