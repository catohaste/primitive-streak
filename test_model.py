import unittest

def test_initial_inhibitor_conditions(delay_steps,inhibitor_low_conc,inhibitor_high_conc,inhibitor):
    for j in range(delay_steps):
        assert(i > inhibitor_low_conc for i in inhibitor[j,:])
        assert(i < inhibitor_high_conc for i in inhibitor[j,:])
