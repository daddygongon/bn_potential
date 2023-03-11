require 'command_line/global'
require 'yaml'

trace = {x: [], y: []}
p vols = [*0..9].collect{|v| 2.5+v/10.0*2.0}
p vols = [*0..9].collect{|v| (2.5+v/10.0*2.0)/5.1279412056 }

vols.each do |val|
  lines = File.readlines("POSCAR_BN")
  lines[1] = "#{val}\n"
  File.write('POSCAR', lines.join(""))
  res = command_line "~/nap/nappy/napsys.py  convert POSCAR pmdini"
  res = command_line "~/nap/pmd/pmd"
  res = command_line "cat erg.pmd"
  p energy = res.stdout.to_f
  trace[:x] << val
  trace[:y] << energy
end

File.write("ewald_res.yaml", YAML.dump(trace))

trace = YAML.load(File.read('ewald_res.yaml'))
trace[:y].map!{|v| v/24.0}
trace[:x].map!{|v| v/3.626002-1.0} # 0.98879
p trace

