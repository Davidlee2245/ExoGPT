"""
Patch to fix DGL import issue with torchdata.datapipes
This should be imported before importing dgl or rfdiffusion
"""
import sys
import types

# Create stub modules for torchdata.datapipes
if 'torchdata.datapipes' not in sys.modules:
    datapipes_module = types.ModuleType('torchdata.datapipes')
    sys.modules['torchdata.datapipes'] = datapipes_module
    
    iter_module = types.ModuleType('torchdata.datapipes.iter')
    sys.modules['torchdata.datapipes.iter'] = iter_module
    
    # Create a stub IterDataPipe class
    class IterDataPipeStub:
        pass
    
    iter_module.IterDataPipe = IterDataPipeStub

