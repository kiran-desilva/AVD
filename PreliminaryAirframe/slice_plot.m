%fs ts hs n w

figure
hold on

fs_idx = [1 2 3 4 5 6 7 8 9];
subplot_dimentions = [3,3];

for i=1:length(fs_idx)
    subplot(subplot_dimentions(1),subplot_dimentions(2),i);
    plot_slice(fs_idx(i),analysis);
    c = colorbar;
    colormap(turbo)
    title(num2str(fs_idx(i)));
    
    % xlabel('N')
    % ylabel('Stringer Thickness')
    % zlabel('Stringer Height')
    % c.Label.String = "Weight";
    % caxis([0 100])
    xlabel('Stringer Thickness')
    ylabel('Stringer Height')
    zlabel('Weight')
    c.Label.String = "N stringers";
    caxis([0 100])
    
    grid on
    
end






function [slice] = create_slice(fs_idx,data)
    flatten = @(x) x(:);

    % slice.x = data.valid_domain.N(fs_idx,:,:,:);
    % slice.y = data.valid_domain.TS(fs_idx,:,:,:);
    % slice.z = data.valid_domain.HS(fs_idx,:,:,:);
    % slice.c = data.valid_domain.weights(fs_idx,:,:,:);
    slice.c = data.valid_domain.N(fs_idx,:,:,:);
    slice.x = data.valid_domain.TS(fs_idx,:,:,:);
    slice.y = data.valid_domain.HS(fs_idx,:,:,:);
    slice.z = data.valid_domain.weights(fs_idx,:,:,:);
end

function [slice,ax] = plot_slice(fs_idx,data)
    slice = create_slice(fs_idx,data);
    ax = scatter3(slice.x(:),slice.y(:),slice.z(:),50,slice.c(:),'filled');
    % ax = surf(slice.x,slice.y,slice.z,slice.c)
end