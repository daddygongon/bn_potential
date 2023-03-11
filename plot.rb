# -*- coding: utf-8 -*-
require 'rbplotly'
require 'yaml'
require './poscar'
require './eam'
require './jindo_lj'

class LJBN < EAM
  attr_accessor :cut_off
  #  <2023-03-10 金>
  # derive LJ in ewald_fitting.mw
AA = 0.932307
BB = 1.480392
VN = 7.761241
VM = 3.880620
R0 = 1.570105

  K2eV = 8.617385e-05
  ALat = 3.626

  def initialize(system, poq='')
    @cut_off = R0*1.3 
#    p @cut_off
    super(system, poq)
  end
  include Math
  def atom_energy(ai)
    rep, rho = 0.0, 0.0
    ai.nl.each do |j|
      r = distance(ai.pos, @system.atoms[j].pos)
      rh = R0/r
      rep += AA * rh ** VN
      rho += BB * rh ** VM
    end
    bind = - rho
    [rep+bind, rep, bind]
  end
  def puts_each_atom_energy()
    @system.atoms[0..9].each_with_index do |atom, i |
      e, _r, _b = atom_energy(atom)
#      p [i, e]
    end
    return ''
  end
end

def e_ewald(x)
  val = 15.73855996*x**3 - 7.697594308*x**2 + 7.500302033*x - 7.511357132 #p8
#  val = 15.73855996*x^3 - 7.697594308*x^2 + 7.500302033*x - 7.511357132 #p24
  val = 76.36350686*x**3 - 61.42859654*x**2 + 60.47310556*x - 66.80819236
  val * 8
end

trace = YAML.load(File.read('ewald_res.yaml'))
trace[:y].map!{|v| v/8.0}
trace[:x].map!{|v| v/3.626002-1.0} # 0.98879
p trace

xx, y1, y2, y3 = [], [], [], []
vols = [*0..10].collect{|v| -0.04 + v*0.01}
vols.each do |val|
  #  vol = 3.626002*val
  x = 1.0 +val 
  xx << x
  ModPoscar.new(x, 0, 0, 0.0)
  poscar = Poscar.new("POSCAR")
  p n_atom = 1 || poscar.n_atoms
  lj = LJBN.new(poscar)
  e_lj = lj.total_energy

  p [val, x, e_ewald(val)/n_atom, e_lj/n_atom, (e_ewald(val) + e_lj)/n_atom]
  y1 << e_ewald(val)
  y2 << e_lj
  y3 << (e_ewald(val) + e_lj)/n_atom
end

traces = [#{x: xx, y: y1},
          #{x: xx, y: y2},
          {x: xx, y: y3}]  
layout = { title: 'Ewald sum',
           xaxis: { title: 'lattice expansion [l/l0]' },
           yaxis: { title: 'energy [eV/system]' } }
pl = Plotly::Plot.new(data: traces,
                      layout: layout)
pl.generate_html(path: './line_chart1.html')


