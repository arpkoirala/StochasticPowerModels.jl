using RepresentativePeriodsFinder, JuMP, Cbc
config_file="C:/Users/karpan/Documents/After_phd/timeslice/default.yaml"
pf = PeriodsFinder(config_file;  populate_entries=true)
optimizer = optimizer_with_attributes(Cbc.Optimizer, "seconds" => 300) # seconds specifies time out limit
pf = find_representative_periods(pf, optimizer=optimizer)