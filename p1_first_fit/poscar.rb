require 'scanf'
require 'matrix'
include Math

#DEBUG = { verbose: true, pos_save: :period_100000 } # pos_save: :p_1 for all, :p_1000000 for never

class Atom
  attr_accessor :pos, :pos0, :nl, :ks, :e0
  def initialize(pos = [0.0, 0.0, 0.0])
    @pos = Vector[0.0,0.0,0.0]
    @pos0 = Vector[0.0,0.0,0.0]
    3.times{|i| @pos0[i] = @pos[i] = pos[i]}
    @nl = []
    @ks = [0.0,0.0,0.0]
    @e0 = 0.0
  end
end

# Modify Poscar for Einstein Calc.
# INPUT: POSCAR_orig
# OUTPU: POSCAR
# with parameters of vol, site, idx, and dev(iation)
class ModPoscar
  def initialize(vol, site=0, idx=0, dev=0.0)
    read_poscar
#    @offset = @lines[7].include?('Direct') ? 8 : 9
    #    p @offset
    volume(vol)
    mod_site(site,idx,dev)
    File.write('POSCAR', @lines.join)
  end

  def read_poscar(source='POSCAR_orig')
    @lines = File.readlines(source)
    @offset = 6
    @dynamics = ' '
    @lines[5..9].each_with_index do |line, i|
      @dynamics = ' T T T' if line.include?('dynamics')
      @offset += i if line.include?('Direct')
    end
  end

  def mod_site(site, idx, dev)
    puts "* fix calc kpoints:50, site:#{site}, xyz_idx:#{idx}, dev:#{dev}"
#    puts "** start: #{DateTime.now}"
    puts calc_a1
    p r_dev = calc_dev(idx, dev)
    p xyz = @lines[site+@offset].scanf("%f %f %f")
    p idx
    xyz[idx] =xyz[idx]+r_dev
    xyz << @dynamics
    @lines[site+@offset] = sprintf("%20.15f %20.15f %20.15f %s\n", *xyz)
  end

  def volume(vol)
    whole = @lines[1].scanf("%f")[0]
    @lines[1] = sprintf("%20.15f\n", whole*vol.to_f)

  end

  def calc_dev(idx, dev)
    whole = @lines[1].scanf("%f")[0]
    l_xyz = @lines[idx+2].scanf("%f %f %f")[idx]
    return dev/l_xyz/whole
  end

  def calc_a1
    whole = @lines[1].scanf("%f")[0]
    l_x = @lines[2].scanf("%f %f %f")[0]
    return "** a1: "+ (whole*l_x/3.0/sqrt(2.0)).to_s
  end
end

class Poscar
  attr_accessor :la, :ions, :element, :title, :atoms, :n_atoms, :lat_vec
  def initialize(file)
    @la, @ions = read_poscar(file)
    @lat_vec = [@la[0][0],@la[1][1],@la[2][2]]
    mk_atoms
  end

  def report(opts)
    case opts[:key]
    when :pos_comp
      @atoms.each do |atom|
        printf("%10.5f %10.5f %10.5f -", *atom.pos0)
        printf("%10.5f %10.5f %10.5f\n", *atom.pos)
      end
    end
  end
  def read_poscar(file)
    lines = File.readlines(file)
    nums = lines[6].scanf("%d %d %d")
    nums = lines[5].scanf("%d %d %d") if nums.size == 0
    @element = lines[7]
    atom = nums.inject(0){|s, i| s+=i }
    la = Array.new(3){Vector[0.0,0.0,0.0]}
    ions = Array.new(atom){Vector[0.0,0.0,0.0]}
    section = []
    index = 0
    lines.each_with_index do |line, i|
      @title = line.chomp if i==0
      if i== 1
        @whole_scale = line.chomp.to_f
        section.push 'la'
        next
      end
      if i==5
        section.pop
        next
      end

      case line
      when /^[D|d]irect/
        section.push 'ions'
        index = 0
        next
      when /^\s*\n/
#        p [i,line] if DEBUG[:verbose]
        break
      end

      case section
      when ["la"]
        line.split(' ').each_with_index{|ele, i| la[index][i] = ele.to_f*@whole_scale }
        index += 1
      when ['ions']
        line.split(' ')[0..2].each_with_index{|ele, i| ions[index][i] = ele.to_f}
        index += 1
      else
      end
    end
    return la, ions
  end
  def expand_la(i, j, k)
    @la[0]*i + @la[1]*j + @la[2]*k
  end
  def expand(nx, ny, nz)
    @la[0][0] = @la[0][0]*nx
    @la[1][1] = @la[1][1]*ny
    @la[2][2] = @la[2][2]*nz
    new_ions = []
    nx.times do |i|
      ny.times do |j|
        nz.times do |k|
          ions.each do |ion|
            new_ion = DFloat::zeros(3)
            new_ion[0] = (i+ion[0])/nx
            new_ion[1] = (j+ion[1])/ny
            new_ion[2] = (k+ion[2])/nz
            new_ions << new_ion
          end
        end
      end
    end
    @ions = new_ions
  end
  def mk_atoms
    @n_atoms = @ions.size
    @atoms = []
    @ions.each do |ion|
      tmp = Vector[0.0,0.0,0.0]
      3.times{|i| tmp[i] = ion[i]*@la[i][i] }
      @atoms << Atom.new(tmp)
    end
  end

  def write_poscar(file, message = "poscar")
    cont = message+"\n"
    cont << "#{@whole_scale}\n"
    @la.each do |lattice|
      cont << sprintf("%15.10f %15.10f %15.10f\n",
                      lattice[0]/@whole_scale,
                      lattice[1]/@whole_scale,
                      lattice[2]/@whole_scale)
    end
    cont << "#{ions.size}\n"
    cont << "Selective dynamics\nDirect\n"
    @atoms.each do |ion|
      cont << sprintf("%15.10f %15.10f %15.10f T T T\n",
                      ion.pos[0]/@la[0][0],
                      ion.pos[1]/@la[1][1],
                      ion.pos[2]/@la[2][2])
    end
    File.write(file, cont)
  end

  def distance(ipos, jpos)
    tmp = 0.0
    3.times do |i|
      x1 = ipos[i] - jpos[i]
      x = x1 - (x1/@lat_vec[i]).round * @lat_vec[i]
      tmp += x * x
    end
    Math.sqrt(tmp)
  end

end

