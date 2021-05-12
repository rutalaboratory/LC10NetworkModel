### LC10a Spiking Network Model

Sample code for simulating the LC10a network model in Hindmarsh Sten et al., 2021. 

LC10a neurons are modeled as a population of integrate-and-fire units, each of which covers an equal portion of the visual field of male flies (estimated as 270° with 15° of binocular overlap). The motion-based receptive fields estimated by Ribeiro et al. (Cell, 2018) are used to structure excitatory input to these model neurons. The total number of spikes of the LC10a population in the right and left hemisphere is integrated and compared by downstream units such that the strength of ipsilateral turning depended on the difference in firing rate between hemispheres in a given time-bin, approximating lateralized input to descending pathways. Importantly, the model is sensitive to the change in the angular position of moving objects, but not their angular size on the retina.

There exists two variants of this model, one tailed for simulating tethered courtship behavior, and the other for simulating the turning behavior of freely courting males given the relative position of the female target. In each case, sample behavioral and stimulus data is provided. 

Requires [Circular Statistics Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
