require_relative './eam'

class  Jindo < EAM
  attr_accessor :cut_off
  #  'Cu' => [
    K2eV = 8.617385e-05
    D0 = 4125.7/2 # * K2eV
    VN = 9.0
    VM = 5.5
    R0 = 2.5487 # e-08
    ALat = 3.6153 # e-08

  def initialize(system, poq)
    @cut_off = 3.6153 * 1.11 #1.07 original # for Jindo's Cu
    super(system, poq)
  end
  include Math
  def atom_energy(ai)
    rep, rho = 0.0, 0.0
    ai.nl.each do |j|
      r = distance(ai.pos, @system.atoms[j].pos)
      rh = R0 / r
      rep += D0 / (VN - VM) * VM * rh ** VN
      rho += D0 / (VN - VM) * VN * rh ** VM
    end
    bind = - rho
    [rep+bind, rep, bind]
  end
end

class JindoEV < EAM
  # D0 is the only one difference to Jindo
  # but I can't find out agood way to use slight different constant.
  # 2022/04/02
  #  'Cu' => [
    K2eV = 8.617385e-05
    D0 = 4125.7/2 * K2eV
    VN = 9.0
    VM = 5.5
    R0 = 2.5487 # e-08
    ALat = 3.6153 # e-08

  def initialize(system, poq)
    @cut_off = 3.6153 * 1.11 #1.07 original # for Jindo's Cu
    super(system, poq)
  end
  include Math
  def atom_energy(ai)
    rep, rho = 0.0, 0.0
    ai.nl.each do |j|
      r = distance(ai.pos, @system.atoms[j].pos)
      rh = R0 / r
      rep += D0 / (VN - VM) * VM * rh ** VN
      rho += D0 / (VN - VM) * VN * rh ** VM
    end
    bind = - rho
    [rep+bind, rep, bind]
  end
  def puts_each_atom_energy()
    @system.atoms[0..9].each_with_index do |atom, i |
      e, _r, _b = atom_energy(atom)
      p [i, e]
    end
    return ''
  end
end

class LJAl < EAM
  attr_accessor :cut_off
  # <2022/06/20 mon>
  # derive LJ in lj_fitting.mw
  AA = 0.066127
  BB = 0.299374
  VN = 10.458180
  VM = 1.853357
  
  R0 = 2.857701
  K2eV = 8.617385e-05
  ALat = 4.0414

  def initialize(system, poq='')
    @cut_off = ALat * 1.11 #1.07 original # for Jindo's Cu
    p @cut_off
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
      p [i, e]
    end
    return ''
  end
end
