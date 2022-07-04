function output_df = createDF(input_df, name)




if strcmp(name, 'trunk')
    Lateroflexion = input_df(:,1);
    Rotation = input_df(:,2);
    Flexion = input_df(:,3);
    output_df    = [Lateroflexion, Rotation, Flexion];
elseif strcmp(name, 'scapula')
    lateralrotation = input_df(:,1);
    Protraction = input_df(:,2);
    Anteriortilt    = input_df(:,3);

    output_df.(name)    = table(lateralrotation, Protraction, Anteriortilt);

else
    Abduction = input_df(:,1);
    Rotation = input_df(:,2);
    Flexion = input_df(:,3);


output_df.(name)    = table(Abduction, Rotation, Flexion);

end

