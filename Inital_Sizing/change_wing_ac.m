function [] = change_wing_ac(xac)
    load('locations')
    load('tailplane')
    load('wing')
    locations.x_ac_w = xac
    locations.x_wing = xac - wing.Xac_from_tip
    tailplane.horizontal.l = locations.x_ac_h - locations.x_ac_w;
    save('tailplane','tailplane')
    save('wing','wing')
    save('locations','locations')
end
