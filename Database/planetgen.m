%create database structure:
planet=struct;
planet(1).name='Earth';
planet(1).solconst=1367;
planet(1).albedo=0.36;
planet(1).IR=234;
planet(1).radius=6371;
planet(1).grav=398600;

planet(2).name='Moon';
planet(2).solconst=1367;
planet(2).albedo=0.073;
planet(2).IR=430;
planet(2).radius=1737;
planet(2).grav=4904.8695;

planet(3).name='Mars';
planet(3).solconst=589;
planet(3).albedo=0.29;
planet(3).IR=390;
planet(3).radius=3396.2;
planet(3).grav=42828.37;

save 'planet.mat' planet

clear
load planet.mat planet
