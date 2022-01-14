clear all
% close all
% home

d = 330e-9; 
r = d/2;

Cu_sint_deg = asind(1/(2*4));

%Cu_volume_full = ((4/6)*pi*(r.^3))
%Cu_cut_off = (pi/3)*((r-(r*cosd(Cu_sint_deg)))^2)*(3*r-(r-(r*cosd(Cu_sint_deg))))
%rat=Cu_cut_off/Cu_volume_full;
Cu_volume = ((4/6)*pi*(r.^3)) - ( (pi/3)*((r-(r*cosd(Cu_sint_deg)))^2)*(3*r-(r-(r*cosd(Cu_sint_deg)))) );
Cu_contact_area =  pi*(r*sind(Cu_sint_deg)).^2;

Cu_resistance = 1.68e-8;
Ag_resistance = 1.59e-8;

Ag_sint_deg = asind(1/(2*2));

no_of_points = 100;
no_of_steps = no_of_points-1;
interval = (Ag_sint_deg-Cu_sint_deg)/no_of_steps;

theta_for_plotting =  zeros(no_of_points,1,1);

Ag_contact_area = zeros(no_of_points,1,1);
final_segment_Ag_faction = zeros(no_of_points,1,1);
segment_final_Ag_volume = zeros(no_of_points,1,1);
segment_final_volume = zeros(no_of_points,1,1);
segment_final_carea = zeros(no_of_points,1,1);
segment_final_length = zeros(no_of_points,1,1);


Ag_Cu_mass_ratio = zeros(no_of_points,1,1);
contact_area_increase = zeros(no_of_points,1,1);

segment_Cu_resistance = zeros(no_of_points,1,1);
final_segment_resistance = zeros(no_of_points,1,1);
resistance_overall = zeros(no_of_points,1,1);
segment_Cu_average_carea = zeros(no_of_points,1,1);
segment_Cu_length = zeros(no_of_points,1,1);

res_av = zeros(no_of_points,1,1);
res_improvement = zeros(no_of_points,1,1);

final_segment_resistance_new = zeros(no_of_points,1,1);
final_segment_resistance2 = zeros(no_of_points,1,1);
resistance_overall2 = zeros(no_of_points,1,1);

noAg_segment_Cu_length = r*cosd(Cu_sint_deg);

fun = @(L) L/((r*sind(acosd(L/r))).^2);
    xmin = 0;
    xmax = (r*cosd(Cu_sint_deg));
    int_allCu_noAg = integral(fun,xmin,xmax,'ArrayValued',true); 
    
    noAg_segment_Cu_resistance = Cu_resistance * 2*int_allCu_noAg / pi;
    
    

int_Cu = zeros(no_of_points,1,1);
int_f_Cu = zeros(no_of_points,1,1);
int_f_Ag = zeros(no_of_points,1,1);

counter = 0;

for theta = Cu_sint_deg:interval:Ag_sint_deg % Cu_sint_deg
    counter = counter+1;
    theta_for_plotting(counter) = theta;
    
    segment_Cu_length(counter) = r*cosd(theta);
    
    segment_final_length(counter) = (r*cosd(Cu_sint_deg)) - (r*cosd(theta));
    segment_final_carea(counter) = pi*(r*sind(theta)).^2;
    segment_final_Cu_volume = (pi/6) * segment_final_length(counter) * ( 3*(r*sind(theta)).^2 + 3*(r*sind(Cu_sint_deg)).^2 + segment_final_length(counter).^2 );
    segment_final_volume(counter) = segment_final_length(counter) * segment_final_carea(counter);
    segment_final_Ag_volume(counter) = segment_final_volume(counter) - segment_final_Cu_volume;
    
    Ag_contact_area(counter) = segment_final_carea(counter)-Cu_contact_area;
    contact_area_increase(counter) = (segment_final_carea(counter)/Cu_contact_area)-1;
        
    Ag_Cu_mass_ratio(counter) = (10490*segment_final_Ag_volume(counter))/(8960*Cu_volume);
            
    if theta ~= Cu_sint_deg % ensure it won't be "NaN"
        final_segment_Ag_faction(counter) = segment_final_Ag_volume(counter) / segment_final_volume(counter) ;
    else
        final_segment_Ag_faction(counter) = 0;
    end
    Cu_faction_in_final_segment = 1-final_segment_Ag_faction(counter); 
    
    fun = @(L) L/((r*sind(acosd(L/r))).^2);
    xmin = 0;
    xmax = segment_Cu_length(counter);
    int_Cu(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    
    segment_Cu_resistance(counter) = Cu_resistance * 2*int_Cu(counter) / pi;
    
    fun = @(L) L/((r*sind(acosd((r*cosd(Cu_sint_deg)-L)/r))).^2);
    xmin = segment_Cu_length(counter);
    xmax = (r*cosd(Cu_sint_deg));
    int_f_Cu(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    
    resistance_segment_final_Cu_part = Cu_resistance * 2*int_f_Cu(counter) / pi;
    
    fun = @(L) L/(((r*sind(acosd((r*cosd(Cu_sint_deg)-L)/r))).^2)-((r*sind(Cu_sint_deg)).^2));
    xmin = segment_Cu_length(counter);
    xmax = (r*cosd(Cu_sint_deg));
    int_f_Ag(counter) = integral(fun,xmin,xmax,'ArrayValued',true); 
    
    resistance_segment_final_Ag_part = Ag_resistance * 2*int_f_Ag(counter) / pi;
    
    final_segment_resistance(counter) = ( (1/resistance_segment_final_Cu_part) + (1/resistance_segment_final_Ag_part) )^(-1);
    resistance_overall(counter) = segment_Cu_resistance(counter) + final_segment_resistance(counter);
   
    res_improvement(counter) = resistance_overall(counter)/noAg_segment_Cu_resistance;
end

%figure,plot(theta_for_plotting,resistance_overall)
%figure,plot(Ag_Cu_mass_ratio,resistance_overall)
figure(1)
hold on
plot(Ag_Cu_mass_ratio,res_improvement)
xlabel('Ag/Cu mass ratio')
ylabel('Resistance reduction due to Ag accumulation in Cu neck')

% figure(2)
% hold on
% plot(Ag_Cu_mass_ratio,resistance_overall)
% xlabel('Ag contact angle (degC)')
% ylabel('Resistance (\Omega)')
