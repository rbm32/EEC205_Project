function r = wallshot_distance(r1,r2, wall_dist)
    %calculates radial distance between two points r1 and r2 after
    %reflection from Z = bar_dist

    % find the intersection point where r1, r2, and wall_dist
    % in the wallshot, the array is assumed to be pointing up such that
    % the distance from antennas to the wall is wall_dist in the +Z
    % direction
    
    % from angle symmetry, the intersection point will be the average of
    % X and Y and wall_dist
    assert(wall_dist > 0, "wall_dist must be greater than 0.");
    intersection_point = (r1 + r2) ./ 2;
    r1(3)=0;
    r2(3)=0;
    intersection_point(3) = wall_dist;
    
    r = radial_distance(r1,intersection_point) ...
        + radial_distance(intersection_point, r2);
end

