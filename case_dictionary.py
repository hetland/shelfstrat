
import numpy as np
import os



class case_dictionary(object):
    
    def __init__(self, directory):
        cases = os.walk(directory).next()[1]
        
        case_dict = {}
        
        for rootdir in cases:
            case = rootdir.split('_')[1:]  # remove leading description 'shelfstrat'
            keys = case[::2]
            vals = [float(val) for val in case[1::2]]
            params = dict(zip(keys, vals))
            params['Ri'] = params['N2'] * params['f']**2 / params['M2']**2
            params['S'] = params['M2'] * 1e-3 * np.sqrt(params['Ri']) / params['f']**2
            params['delta'] = np.sqrt(params['Ri']) * params['S']
            
            case_dict[rootdir] = params
            
        self._case_dict = case_dict
    
    def find(self, **args):
        
        cases = []
        
        for case, case_params in self._case_dict.iteritems():
            
            matches = {p: False for p in args.keys()}
            
            for case_param, case_param_val in case_params.iteritems():
                for trial_param, trial_val in args.iteritems():
                    if case_param == trial_param and np.isclose(case_param_val, trial_val, rtol=0.01):
                        matches[trial_param] = True
            
            if np.alltrue(matches.values()):
                cases.append(case)
        
        return cases
    
    def values(self, param):
        values = []
        for case, case_params in self._case_dict.iteritems():
            values.append(case_params[param])
        
        return values
    
    def __getitem__(self, case):
        return self._case_dict[case]


if __name__ == '__main__':
    cases = case_dictionary('./simulations')
    matches = cases.find(delta=0.3)
    for match in matches:
        print match, '  ', cases[match]