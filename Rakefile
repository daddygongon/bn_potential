# -*- coding: utf-8 -*-
require "colorize"
require 'command_line/global'

task :default do
  system "rake -T"
end

desc 'git auto'
task :git_auto do
  print "Input comments: "
  comment = STDIN.gets.chomp 
  comment = comment == "" ? "auto pull and push" : comment.gsub("\'", "\\'")
  ["git add -A",
   "git commit -m '#{comment}'",
   "git pull origin main",
   "git push origin main"].each do |comm|
    puts comm.cyan
    res = command_line comm
    puts res.stderr.red
    puts res.stdout.green
  end
end

