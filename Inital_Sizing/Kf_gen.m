Kf_data = [8.895423400387589, 0.05462753950338595
18.89282646477168, 0.25477802859292664
24.849078087856473, 0.4675696012039121
30.702011626281166, 0.6740406320541754
60.65514692662693, 2.9936794582392774
43.60095087796425, 1.0322046651617751];

Kf_fit_data = [61.65860284863898, 3.122197140707298
60.47248247068495, 2.9283671933784796
58.47324157494191, 2.6544770504138446
57.26410835214482, 2.5112114371708047
56.27311772108945, 2.35530474040632
55.37010327813212, 2.2394281414597437
54.27004135120589, 2.0898419864559816
53.05036057452237, 1.969751693002257
52.23436345112799, 1.8960120391271629
50.600451467268925, 1.752746425884123
48.77716294772199, 1.5589164785553042
47.12982680437112, 1.445146726862302
45.89959847380119, 1.3482317531978927
44.248426856309806, 1.2428893905191867
42.79813819692788, 1.1628291948833707
41.556403443936404, 1.0911963882618507
39.07101620088319, 0.9521444695259587
37.19594878043931, 0.8720842738901422
34.695219640824035, 0.7667419112114366
33.44389620248127, 0.716177577125658
31.356439401506403, 0.6361173814898415
28.633732196008882, 0.5518434913468768
26.75099382728396, 0.4886380737396534
24.655866078028073, 0.42543265613242953
22.560738328772178, 0.3622272385252061
21.091272298687674, 0.3243039879608718
17.41568948640614, 0.2337095560571858
13.946742843444813, 0.15575620767494325
11.000139834994762, 0.09676448457486808
9.000659222117946, 0.056734386756959854
];

Kf_fit = fit(Kf_fit_data(:,1),Kf_fit_data(:,2),'poly4');

figure
hold on
plot(Kf_data(:,1),Kf_data(:,2),'x');
plot(Kf_fit_data(:,1),Kf_fit_data(:,2),'o')
plot(Kf_fit);
xlabel('wing quater chord position, percent of fuselage length')
ylabel('Kf')