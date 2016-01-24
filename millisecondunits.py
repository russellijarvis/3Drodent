grep "msecond" units.py
    "msecond",
    "msecond2",
    "msecond3",
msecond = Unit.create_scaled_unit(second, "m")
msecond2 = msecond**2
msecond2.name = "msecond2"
msecond3 = msecond**3
msecond3.name = "msecond3"
    msecond = float(msecond)
    msecond2 = float(msecond2)
    msecond3 = float(msecond3)
    msecond,
    msecond2,
    msecond3,
    msecond,
    msecond2,
    msecond3,
    pE("", "get_unit(3*msecond)")

