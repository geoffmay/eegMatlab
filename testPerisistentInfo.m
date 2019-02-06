function testPerisistentInfo()

defaultArgs = {10, 'trainlm'};

load('C:\Users\Neuro\Documents\MATLAB\resources\defaultNeuralNetParameters.mat');
[args,param] = nnparam.extract_param(defaultArgs,INFO.defaultParam);
[param,err] = INFO.overrideStructure(param,args);
if ~isempty(err), nnerr.throw('Args',err,'Parameters'); end
net = create_network(param);
net.name = INFO.name;
out1 = init(net);


end