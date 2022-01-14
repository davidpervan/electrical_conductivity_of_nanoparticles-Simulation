clear all
% close all
% home

l_d = 330e-9; 
l_r = l_d/2;
particle_size_ratio = 1/1.0;
r_d = l_d*particle_size_ratio;
r_r = r_d/2;

l_Cu_sint_deg = asind(1/(2*8)); % must be smaller than r_r

if l_r*sind(l_Cu_sint_deg) >= r_r
    l_Cu_sint_deg = asind((r_r*.999)/l_r)
end

r_Cu_sint_deg = asind(l_r*sind(l_Cu_sint_deg)/r_r);

%Cu_volume_full = ((4/6)*pi*(r.^3))
%Cu_cut_off = (pi/3)*((r-(r*cosd(Cu_sint_deg)))^2)*(3*r-(r-(r*cosd(Cu_sint_deg))))
%rat=Cu_cut_off/Cu_volume_full;
l_Cu_volume = ((4/6)*pi*(l_r.^3)) - ( (pi/3)*((l_r-(l_r*cosd(l_Cu_sint_deg)))^2)*(3*l_r-(l_r-(l_r*cosd(l_Cu_sint_deg)))) );
r_Cu_volume = ((4/6)*pi*(r_r.^3)) - ( (pi/3)*((r_r-(r_r*cosd(r_Cu_sint_deg)))^2)*(3*r_r-(r_r-(r_r*cosd(r_Cu_sint_deg)))) );
both_Cu_contact_area =  pi*(l_r*sind(l_Cu_sint_deg)).^2;

Cu_resistance = 1.68e-8;
Ag_resistance = 1.59e-8;

l_Ag_sint_deg = asind(1/(2*3)); % must be smaller than r_r hence must be smaller or equal to particle size ratio % must be bigger than l_Cu_sint_deg

if l_r*sind(l_Ag_sint_deg) > r_r
    l_Ag_sint_deg = asind(r_r/l_r)
end

r_Ag_sint_deg = asind(l_r*sind(l_Ag_sint_deg)/r_r);

no_of_points = 100;
no_of_steps = no_of_points-1;
interval = (l_Ag_sint_deg-l_Cu_sint_deg)/no_of_steps;

l_theta_for_plotting =  zeros(no_of_points,1,1);
r_theta_for_plotting =  zeros(no_of_points,1,1);

Ag_contact_area = zeros(no_of_points,1,1);
l_final_segment_Ag_faction = zeros(no_of_points,1,1);
r_final_segment_Ag_faction = zeros(no_of_points,1,1);

l_segment_final_Ag_volume = zeros(no_of_points,1,1);
l_segment_final_volume = zeros(no_of_points,1,1);
both_segment_final_carea = zeros(no_of_points,1,1);
l_segment_final_length = zeros(no_of_points,1,1);

r_segment_final_length = zeros(no_of_points,1,1);
r_segment_final_volume = zeros(no_of_points,1,1);
r_segment_final_Ag_volume = zeros(no_of_points,1,1);

Ag_Cu_mass_ratio = zeros(no_of_points,1,1);
contact_area_increase = zeros(no_of_points,1,1);

l_segment_Cu_resistance = zeros(no_of_points,1,1);
r_segment_Cu_resistance = zeros(no_of_points,1,1);
l_final_segment_resistance = zeros(no_of_points,1,1);
r_final_segment_resistance = zeros(no_of_points,1,1);

resistance_overall = zeros(no_of_points,1,1);
segment_Cu_average_carea = zeros(no_of_points,1,1);
l_segment_Cu_length = zeros(no_of_points,1,1);
r_segment_Cu_length = zeros(no_of_points,1,1);

res_av = zeros(no_of_points,1,1);
res_improvement = zeros(no_of_points,1,1);

final_segment_resistance_new = zeros(no_of_points,1,1);
final_segment_resistance2 = zeros(no_of_points,1,1);
resistance_overall2 = zeros(no_of_points,1,1);

fun = @(L) L/((l_r*sind(acosd(L/l_r))).^2);
    xmin = 0;
    xmax = (l_r*cosd(l_Cu_sint_deg));
    l_int_allCu_noAg = integral(fun,xmin,xmax,'ArrayValued',true); 
    l_noAg_segment_Cu_resistance = Cu_resistance * l_int_allCu_noAg / pi;

fun = @(L) L/((r_r*sind(acosd(L/r_r))).^2);
    xmin = 0;
    xmax = (r_r*cosd(r_Cu_sint_deg));
    r_int_allCu_noAg = integral(fun,xmin,xmax,'ArrayValued',true); 
    r_noAg_segment_Cu_resistance = Cu_resistance * r_int_allCu_noAg / pi;
    
overall_noAg_segment_Cu_resistance = l_noAg_segment_Cu_resistance + r_noAg_segment_Cu_resistance;

l_int_Cu = zeros(no_of_points,1,1);
l_int_f_Cu = zeros(no_of_points,1,1);
l_int_f_Ag = zeros(no_of_points,1,1);

r_int_Cu = zeros(no_of_points,1,1);
r_int_f_Cu = zeros(no_of_points,1,1);
r_int_f_Ag = zeros(no_of_points,1,1);

counter = 0;

for l_theta = l_Cu_sint_deg:interval:l_Ag_sint_deg % Cu_sint_deg
    r_theta = asind(l_r*sind(l_theta)/r_r);
    counter = counter+1;
    l_theta_for_plotting(counter) = l_theta;
    r_theta_for_plotting(counter) = r_theta;
    
    l_segment_Cu_length(counter) = l_r*cosd(l_theta);
    r_segment_Cu_length(counter) = r_r*cosd(r_theta);
    
    l_segment_final_length(counter) = (l_r*cosd(l_Cu_sint_deg)) - (l_r*cosd(l_theta));
    both_segment_final_carea(counter) = pi*(l_r*sind(l_theta)).^2;
    l_segment_final_Cu_volume = (pi/6) * l_segment_final_length(counter) * ( 3*(l_r*sind(l_theta)).^2 + 3*(l_r*sind(l_Cu_sint_deg)).^2 + l_segment_final_length(counter).^2 );
    l_segment_final_volume(counter) = l_segment_final_length(counter) * both_segment_final_carea(counter);
    l_segment_final_Ag_volume(counter) = l_segment_final_volume(counter) - l_segment_final_Cu_volume;
    
    r_segment_final_length(counter) = (r_r*cosd(r_Cu_sint_deg)) - (r_r*cosd(r_theta));
    r_segment_final_Cu_volume = (pi/6) * r_segment_final_length(counter) * ( 3*(r_r*sind(r_theta)).^2 + 3*(r_r*sind(r_Cu_sint_deg)).^2 + r_segment_final_length(counter).^2 );
    r_segment_final_volume(counter) = r_segment_final_length(counter) * both_segment_final_carea(counter);
    r_segment_final_Ag_volume(counter) = r_segment_final_volume(counter) - r_segment_final_Cu_volume;
        
    Ag_contact_area(counter) = both_segment_final_carea(counter)-both_Cu_contact_area;
    contact_area_increase(counter) = (both_segment_final_carea(counter)/both_Cu_contact_area)-1;
        
    Ag_Cu_mass_ratio(counter) = (10490*(l_segment_final_Ag_volume(counter)+r_segment_final_Ag_volume(counter)))/(8960*(l_Cu_volume+r_Cu_volume));
            
    if l_theta ~= l_Cu_sint_deg % ensure it won't be "NaN"
        l_final_segment_Ag_faction(counter) = l_segment_final_Ag_volume(counter) / l_segment_final_volume(counter) ;
    else
        l_final_segment_Ag_faction(counter) = 0;
    end
    l_Cu_faction_in_final_segment = 1-l_final_segment_Ag_faction(counter); 
    
    if r_theta ~= r_Cu_sint_deg % ensure it won't be "NaN"
        r_final_segment_Ag_faction(counter) = r_segment_final_Ag_volume(counter) / r_segment_final_volume(counter) ;
    else
        r_final_segment_Ag_faction(counter) = 0;
    end
    r_Cu_faction_in_final_segment = 1-r_final_segment_Ag_faction(counter); 
    
    fun = @(L) L/((l_r*sind(acosd(L/l_r))).^2);
    xmin = 0;
    xmax = l_segment_Cu_length(counter);
    l_int_Cu(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    l_segment_Cu_resistance(counter) = Cu_resistance * l_int_Cu(counter) / pi;
    
    fun = @(L) L/((r_r*sind(acosd(L/r_r))).^2);
    xmin = 0;
    xmax = r_segment_Cu_length(counter);
    r_int_Cu(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    r_segment_Cu_resistance(counter) = Cu_resistance * r_int_Cu(counter) / pi;
    
    fun = @(L) L/((l_r*sind(acosd((l_r*cosd(l_Cu_sint_deg)-L)/l_r))).^2);
    xmin = l_segment_Cu_length(counter);
    xmax = (l_r*cosd(l_Cu_sint_deg));
    l_int_f_Cu(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    l_resistance_segment_final_Cu_part = Cu_resistance * l_int_f_Cu(counter) / pi;
    
    fun = @(L) L/((r_r*sind(acosd((r_r*cosd(r_Cu_sint_deg)-L)/r_r))).^2);
    xmin = r_segment_Cu_length(counter);
    xmax = (r_r*cosd(r_Cu_sint_deg));
    r_int_f_Cu(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    r_resistance_segment_final_Cu_part = Cu_resistance * r_int_f_Cu(counter) / pi;
    
    fun = @(L) L/(((l_r*sind(acosd((l_r*cosd(l_Cu_sint_deg)-L)/l_r))).^2)-((l_r*sind(l_Cu_sint_deg)).^2));
    xmin = l_segment_Cu_length(counter);
    xmax = (l_r*cosd(l_Cu_sint_deg));
    l_int_f_Ag(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    l_resistance_segment_final_Ag_part = Ag_resistance * l_int_f_Ag(counter) / pi;
    
    fun = @(L) L/(((r_r*sind(acosd((r_r*cosd(r_Cu_sint_deg)-L)/r_r))).^2)-((r_r*sind(r_Cu_sint_deg)).^2));
    xmin = r_segment_Cu_length(counter);
    xmax = (r_r*cosd(r_Cu_sint_deg));
    r_int_f_Ag(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    r_resistance_segment_final_Ag_part = Ag_resistance * r_int_f_Ag(counter) / pi;
    
    l_final_segment_resistance(counter) = ( (1/l_resistance_segment_final_Cu_part) + (1/l_resistance_segment_final_Ag_part) )^(-1);
    r_final_segment_resistance(counter) = ( (1/r_resistance_segment_final_Cu_part) + (1/r_resistance_segment_final_Ag_part) )^(-1);
    
    resistance_overall(counter) = l_segment_Cu_resistance(counter) + l_final_segment_resistance(counter) + r_segment_Cu_resistance(counter) + r_final_segment_resistance(counter);
   
    res_improvement(counter) = resistance_overall(counter)/overall_noAg_segment_Cu_resistance;
end

%figure,plot(theta_for_plotting,resistance_overall)
%figure,plot(Ag_Cu_mass_ratio,resistance_overall)
figure(1)
plot(Ag_Cu_mass_ratio,res_improvement)
xlabel('Ag/Cu mass ratio')
ylabel('Resistance reduction due to Ag accumulation in Cu neck')
hold on

% figure(2)
% hold on
% plot(Ag_Cu_mass_ratio,resistance_overall)
% xlabel('Ag contact angle (degC)')
% ylabel('Resistance (\Omega)')
