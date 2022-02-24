
function [analysis] = fuselage_search_analysis(data,doplot)
    
    %plot full domian
    analysis.full_domain = get_domain(data);
    % scatter_plot(analysis.full_domain);
    

    %filter for failure codes
    analysis.fail_domain = get_filtered_domian(data,data.weights(:)<0);
    

    %filter for valid solutions
    analysis.valid_domain = get_filtered_domian(data,data.weights(:)>0);
    
    %find the best solution
    [~,analysis.valid_domain.best_idx] = min(analysis.valid_domain.weights(:));
    analysis.best_fuselage = analysis.valid_domain.fuselages{analysis.valid_domain.best_idx};
    analysis.valid_domain.valid_fuselages = cellfun(@(x) x, analysis.valid_domain.fuselages(analysis.valid_domain.idx));
    
    if doplot
        multi_variable_plot(analysis.full_domain);
        multi_variable_plot(analysis.fail_domain);
        multi_variable_plot(analysis.valid_domain);
        %3d slice plots
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
            % xlabel('Fuselage Thickness')
            % ylabel('Stringer Thickness')
            % zlabel('Stringer Height')
            % c.Label.String = "Weight";
            % caxis([200 2000])

            xlabel('Stringer Thickness')
            ylabel('Stringer Height')
            zlabel('Weight')
            c.Label.String = "N stringers";
            caxis([0 100])

            grid on

        end
    end


    function [fig] = multi_variable_plot(data)
        figure
        hold on

        X_new = [data.N(:),data.TS(:),data.HS(:),data.weights(:)];
        varNames = {'N stringers'; 'Stringer Thickness'; 'Stringer Height'; 'weights'};

        gplotmatrix(X_new,[],data.FS(:),[],[],15,[],[],varNames)
    end

    function [domain] = get_domain(data)
        domain.FS = data.FS;
        domain.TS = data.TS;
        domain.HS = data.HS;
        domain.N = data.N;
        domain.weights = data.weights;
        domain.fuselages = data.fuselages;
    end

    function [filtered_data] = get_filtered_domian(data,idx)
        %filter for failure codes
        filtered_data = get_domain(data);
        filtered_data.idx = idx;
        nanidx = ~idx;
        filtered_data.FS(nanidx) = NaN;
        filtered_data.TS(nanidx) = NaN;
        filtered_data.HS(nanidx) = NaN;
        filtered_data.N(nanidx) = NaN;
        filtered_data.weights(nanidx) = NaN;
    end

    function [fig] = scatter_plot(data)
        figure
        hold on
        scatter3(data.TS(:),data.HS(:),data.N(:),100000*data.FS(:),data.weights(:));
        c = colorbar;
        zlabel('N Stringers')
        xlabel('Stringer Thickness [m]')
        ylabel('Stringer Height [m]')
        c.Label.String = ' Total Weight [Kg]'
    end

    function [fig] = multi_variable_failure_plot(data)
        figure
        hold on

        X_new = [data.N(:),data.TS(:),data.HS(:),data.FS(:)];
        varNames = {'N stringers'; 'Stringer Thickness'; 'Stringer Height'; 'Fuselage Thickness'};

        gplotmatrix(X_new,[],data.weights(:),[],[],15,[],[],varNames)
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

end