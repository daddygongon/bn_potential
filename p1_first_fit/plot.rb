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
AA = 0.924754
BB = 1.473289
VN = 7.822700
VM = 3.911350
R0 = 1.552398

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
  x = x/3.626002-1.0
  val = 10.65283172*x**3 + (-1)*7.818569264*x**2 + 7.430530135*x - 7.509501261
  val * 8
end

trace = YAML.load(File.read('ewald_res.yaml'))
trace[:y].map!{|v| v/8.0}
trace[:x].map!{|v| v/3.626002-1.0} # 0.98879
p trace

xx, y1, y2, y3 = [], [], [], []
#vols = [*0..9].collect{|v| -0.04+v*0.005}
vols = [*0..9].collect{|v| -0.2+v*0.05}
vols.each do |val|
  x = 1+val
  xx << x
  ModPoscar.new(x, 0, 0, 0.0)
  poscar = Poscar.new("POSCAR")
  l0 = poscar.lat_vec[0]
  lj = LJBN.new(poscar)
  e_lj = lj.total_energy
  
  y1 << e_ewald(l0)/8
  y2 << e_lj/8
  y3 << (e_ewald(l0) + e_lj)/8
end

traces = [{x: xx, y: y1, name: 'Ewald sum'},
          {x: xx, y: y2, name: 'Lennard-Jones'},
          {x: xx, y: y3, name: 'Total energy'}]  
layout = { title: 'Ewald sum',
           xaxis: { title: 'lattice expansion [l/l0]' },
           yaxis: { title: 'energy [eV/system]' } }
pl = Plotly::Plot.new(data: traces,
                      layout: layout)
pl.generate_html(path: './line_chart1.html')


