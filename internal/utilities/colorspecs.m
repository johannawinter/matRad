function color = colorspecs()

color.mpg = [0,0.4717,0.4604]; % color [0,125,122]
color.dre = [0.4906,0,0]; % color [130,0,0]
color.ora = [255,153,51] ./ 255;
color.blu = [0,0,0.509];
color.gra = 0.5 * ones(3,1);
color.red = [233/255,54/255,22/255];
color.gre = [1/255,126/255,51/255];

color.dkfzlB = [169/255,188/255,213/255]; % color [0,125,122]
color.dkfzmB = [70/255,118/255,173/255]; % color [0,125,122]
color.dkfzdB = [0, 75/255, 142/255]; % color [0,125,122]


color.lightmpg = [1,1,1] - 0.5 * ([1,1,1] - color.mpg);
color.lightdre = [1,1,1] - 0.5 * ([1,1,1] - color.dre);
color.lightblu = [1,1,1] - 0.5 * ([1,1,1] - color.blu);
color.lightora = [1,1,1] - 0.5 * ([1,1,1] - color.ora);

end