# These scripts work for the version of grbllh at
# http://code.icecube.wisc.edu/svn/sandbox/richman/grbllh/trunk -r131724
## a new official project version will be released in summer 2015, with minor user-end changes to the code below


# You have made some sort of event selection on your signal and background events and now would like to test its sensitivity against an ensemble of known GRBs with some assumed signal flux


# background -> I use data events not within 2 hours of any GRB time window (given on GRBweb)
# signal -> I use NuE nugen for BDT training and signal energy PDF; I use NuE, NuTau, and NuMu nugen for signal injection into the likelihood trials for sensitivity calculations (I'm assuming you are using all three, since I imagine angular resolution is comparably bad for each flavor at low energies; if not, just remove what you don't need)


# You can use (1) diffuse E^-2-weighted nugen as signal or (2) pseudo-point-source simulation generated from the diffuse simulation.  In (2), events in an x^degree circle around the zenith,azimuth direction of each GRB are used for that GRB.  The latter is the way to go for assuming any burst-specific model neutrino spectrum for signal.
## Below, I am using (2) for signal injection and (1) for the signal energy PDF.  


# 1. llh_psarrays() to create the pseudo-point-source simulation from diffuse nugen, described in (2) above
# 2. llh_sources_events() to set up grbllh Sources and Events objects for the likelihood analysis
# 3. llh_pdfs() to set up the background space, and signal/background energy probability distribution function grbllh objects
# 4. llh_throwers() to construct the grbllh throwers (that 'thrower' signal and background events from an ensemble of grbs)
# 5. llh_null_tsd() calculate and save the test statistic distribution from the background thrower only
# 6. llh_sensitivity_ps() to calculate the limit setting and discovery potentials for the event selection, using the thrower and nulltsd objects


def llh_psarrays ():
    """Create pseudo-point-source simulation from diffuse nugen"""

    # Wrangle the GRB data:
    ## You can make numpy arrays / dictionaries of the GRB timing and spacial info using the text file output of http://icecube.wisc.edu/~grbweb
    ## Or you can get the data directly from the database on dbs3.icecube.wisc.edu (need to be on a cobalt machine at wisc)
    ## The grblist.py script does this.  It uses this class called Vars that Mike Richman wrote and I also included, which just adds object.key attributes to dictionaries.  It's useful, but you obv. don't need to use it.  The most useful part of the grblist.py script is showing how to read from the grbweb db.
    ## grbs below is an object made with this script.
    

    # Make the pseudo-point-source data from the diffuse nugen simulation sets:
    ## Can use code from the attached fakeps.py script that Mike Richman wrote to do this quickly
    ## I use the Arrays class that Mike also wrote to organize the pssim datasets.  This just adds a few things to the Vars class mentioned above, like checking to make sure all included arrays are of the same length and the ability to copy and append.  Again, this is not required for grbllh use.
    def get_ps (grbs, diffsim, flavor):
        
        # fraction of sky (I use 1%) to define circle width about GRB coords to throw nu events
        ## this could be set to the respective error circles around each GRB, but a simple constant circle of ~11degrees about each GRB gives sufficient statistics
        fos = .01
        pssim = Arrays ()

        # indices to organize the events by GRB
        idxs = []
        pssim.weight_E2 = []  # you can add other flux models too
        grb_ns = [] # store GRB numbers
        grb_names = [] # store GRB names

        for n, grb in enumerate (grbs):
            idx = fakeps.ps_events (grb, diffsim, frac_of_sphere=fos) \
                        * (1e0 <= diffsim.true_energy) \
                        * (diffsim.true_energy <= 1e9)
            ## (can also include energy bounds)

            if flavor == 'nue':
                nu = 66
                nubar = 67

            if flavor == 'numu':
                nu = 14
                nubar = -14

            if flavor == 'nutau':
                nu = 133
                nubar = 134

            for type in (nu, nubar):
                type_idx = (diffsim.type == type)
                keep_idx = type_idx * idx  # nu or nubar events to keep for this GRB
                idxs.append (keep_idx)
                grb_ns.append (n)
                grb_names.append (grb.name)
                # NOTE: 4 pi (ngen / 2) because of half ngen per type for nugen
                pssim.weight_E2 = np.r_[pssim.weight_E2,
                                        diffsim.oneweight[keep_idx]
                                        * diffsim.true_energy[keep_idx]**-2
                                        / (4 * np.pi
                                           * fos * diffsim.meta.n_gen / 2.)]
        for n, k in enumerate (diffsim):
            # 'meta' key in my arrays has ngen (number of nugen events generated)
            if k == 'meta':
                continue
            if k.find ('weight') == 0:
                continue
            # add all nugen events and their variables, ordered by GRB
            ## (events can be and are re-used)
            pssim[k] = np.concatenate (
                tuple (diffsim[k][idx] for idx in idxs))

        pssim.meta.frac_of_sphere = fos
        pssim.meta.diffuse_n_gen = diffsim.meta.n_gen
        pssim.meta.ps_n_gen = fos * diffsim.meta.n_gen  # for record of proper weighting, as done above
        pssim.grb_idx = np.concatenate (
            tuple (grb_n * np.ones (np.sum (idx), dtype=int)
                   for (grb_n, idx) in zip (grb_ns, idxs)))
        pssim.grb_names = np.concatenate (
            tuple (np.array (np.sum (idx) * [grb_name])
                   for (grb_name, idx) in zip (grb_names, idxs)))
        # if GBM-located GRB (2 component error - see grblist.py and arxiv.org/abs/1201.3099)
        pssim.gbm_pos = grbs.arrays.gbm_pos[pssim.grb_idx]
        pssim.grb_sigma = grbs.arrays.sigma[pssim.grb_idx] # per-burst localization error
        pssim.grb_t_100 = grbs.arrays.t_100[pssim.grb_idx]

        return pssim

    #ps_nue = get_ps (grbs, nugen_nue, 'nue')
    ps_numu = get_ps (grbs, diffsim, 'numu')
    #ps_nutau = get_ps (grbs, nugen_nutau, 'nutau')

    #saving (ps_nue, 'ps_nue.arrays')
    save (ps_numu, 'ps_numu.arrays')
    #saving (ps_nutau, 'ps_nutau.arrays')


def llh_sources_events ():
    """Set up grbllh Events and Soures objects"""

    # The Sources class sets up the fake random GRBs, either generated from the simulation for signal or from the actual GRB sample for background, used for the likelihood analysis pseudo-search trials 
    # The Events class sets up the events thrown from the direction of each Source for the likelihood analysis pseudo-serach trials


    # Background
    bg_sources = grbllh.Sources (
        grbs.arrays.zenith, grbs.arrays.azimuth,
        grbs.arrays.sigma, grbs.arrays.t_100,
        t=grbs.arrays.t_start,
        gbm_pos=grbs.arrays.gbm_pos)
        
    # 'exp' is an Arrays object of data events outside of 2 hr time windows about each GRB T100
    bg_events = grbllh.Events (
        exp.zenith, exp.azimuth,
        exp.error, exp.eproxy,
	livetime=28475585.0, t=exp.time,
        run=exp.run)
    #bg_events.scores = exp.scores  # BDT scores


    # Signal
    # for each flavor do:
    ## Using pssim for grbllh thrower objects
    ps_sources = grbllh.Sources (
        pssim.true_zenith, 
        pssim.true_azimuth,
        pssim.grb_sigma,
        pssim.grb_t_100,
        gbm_pos=pssim.gbm_pos,
        source_index=pssim.grb_idx,
        smear=True) # use grb_sigma to smear events throwing about zen,az
    ps_events = grbllh.Events (
        pssim.reco_zenith,
        pssim.reco_azimuth,
        pssim.reco_err, # 'corrected' wrt energy with pull distribution
        pssim.eproxy,
        weights={'E2': pssim.weight_E2,
#		 'E3':  pssim.weight_E3,
#                 'E25':  pssim.weight_E25,
                 'SubPhoto':  pssim.subphoto_weight,
		 'ChkGRB': pssim.chkgrb_weight,})
    #ps_events.scores = pssim.scores # BDT scores

    ## Using diffsim for grbllh S/B energy PDF ratio object
    sig_events_diff = grbllh.Events (
        diffsim.reco_zenith, diffsim.reco_azimuth,
        diffsim.reco_err, diffsim.eproxy,
        weights={
            'E2': diffsim.weight_E2,})
    #sig_events_diff.scores = diffsim.scores


    save (bg_sources, 'bg.sources')
    save (bg_events, 'bg.events')
    save (ps_sources, 'ps.sources')
    save (ps_events, 'ps.events')
    save (sig_events_diff, 'diff.events')


def llh_pdfs ():
    """Set up grbllh Space and Energy PDFs with diffuse simulation"""

    pdf_space_bg = grbllh.PDFSpaceBg (bg_events)
    save (pdf_space_bg, 'pdf_space_bg.pickle')

    # using data here for background energy PDF (statistics up to 1 PeV)
    ## can also use atmospheric-weighted nugen
    #pdf_ratio_energy = grbllh.PDFRatioEnergy (
    #    sig_events_diff, 'E2',
    #    data_events=None,  ###bg_events
    #	bg_weight_name = "E2",
    #    range=(0.5, 4)) # in log10(GeV) (I use 2, 6 from my bg energy range)
    pdf_ratio_energy = grbllh.FlatPDFRatioEnergy()  ### New Flat Energy PDF
    save (pdf_ratio_energy, 'pdf_ratio_energy.pickle')
    

def llh_throwers ():
    """Set up grbllh signal and background thrower objects"""

    # 'r' is the S/B PDF ratio (here, space*time*energy)
    ## typically, 0.1 is a good threshold for semi-significant events
    ## setting the 'min_r' threshold in the grbllh.Config class helps speed up processing
    ## for unblinding, I set 'min_r' equal to 0 (no threshold - so I see all events)
    min_r = 0.1
    # 'sigma_t_truncate' is how many sigmas should the Gaussian tails of the signal time PDF extend on either side of the flat part during the GRB T100
    ## typically, 4. has been set to include early and late emission
    sigma_t_truncate = 4.
    # 'max_delang' between event and grb, I set to pi since my events have uncertainties that high
    max_delang = np.pi

    config = grbllh.Config (
        min_r=min_r,
        sigma_t_truncate=sigma_t_truncate,
        max_delang=max_delang)


    bg_thrower = grbllh.PseudoBgThrower (
		"2012_LowEn",
                bg_events,
                bg_sources,
                pdf_space_bg,
                pdf_ratio_energy,
                energy_bins = 5,
                zenith_bins = 4,  # per energy bin
                config=config)
    save (bg_thrower,
            'bg_thrower.pickle')

    def save_ps_thrower (flavor, events, sources, pdf_space_bg, pdf_ratio_energy):
        for weight_name in sorted (events.weights):
            weights = events.weights[weight_name]
            thrower  = grbllh.SignalThrower ("2012_LowEn",
                events, weight_name, sources,
                pdf_space_bg, pdf_ratio_energy,
                mu=1.0, # multiplying factor on weights
                config=config)
            save (thrower,  'ps_{0}_{1}.thrower'.format (flavor,weight_name))

    #save_ps_thrower ('nue', ps_nue_events, ps_nue_sources, pdf_space_bg, pdf_ratio_energy)
    save_ps_thrower ('numu', ps_events, ps_sources, pdf_space_bg, pdf_ratio_energy)
    #save_ps_thrower ('nutau', ps_nutau_events, ps_nutau_sources, pdf_space_bg, pdf_ratio_energy)


# To calculate sensitivities per BDT cut score or some type of event selection,
# like shown at http://umdgrb.umd.edu/~hellauer/GRBAnalysis/plots/optimization79/lim_discpot_paper-new-pcuts.png,
# create PDFs and throwers for events at each cut level to test
# then can calculate the test statistic distributions, sensitivity and discovery potentials for each cut as below

# See http://icecube.umd.edu/~rmaunu/docs/grb_ts_internal_note/RyanMaunu_GRB_TS_Internal_Note.pdf for different ways to calculate the test statistic
## Below, I am using the 'classic' 'stacked' method
## Changing to the 'max(T_burst') method is fairly straightforward; to do this, see the function options in grbllh/__init__.py
## I typically do 1e8 trials for the null_tsd calculations on the cluster and 5e4 trials for the limit setting and disovery potential calculations on the cluster

def llh_null_tsd ():
    """calculate and save the test statistic distribution from the background thrower only"""

    seed = time.time () * os.getpid ()  # some random number generator seed
    bg_throwers = [bg_thrower]
    tsd = grbllh.do_trials (n_trials, bg_throwers, seed=seed, llh_type=llh_type)
    tsd.T_vals # force sort
    tsd.seed = seed
    save (tsd, 'tsd.pickle')


def llh_sensitivity_ps (beta,thresh,median=False):
    """Calculate limit setting or discovery potential for a given event selection and flux model"""

    sig_throwers = [sig_thrower_numu]
    bg_throwers = [bg_thrower]
    expectation = np.sum ([thrower.prob.sum () for thrower in sig_throwers])
    # prob = event weight (E^-2)
    mu_top = 2. / expectation
    # so mu_top is a little less than Poisson limit (2.3 / expected number of events)

    # Calculating signal strength for T > {0} in pct% of trials
    ## Examples:
    # T_thrsh = 0
    if median:
    	T_thresh = null_tsd.prob_thresh (0.5) ## median value
    else:
	T_thresh = null_tsd.sigma_thresh (thresh)  ## 3, 4, 5 sigma values from null tsd

    pct = beta
    #pct = 0.9
    stuff = Vars (grbllh.calc_mu_beta_thresh (
            T_thresh, pct, n_trials, bg_throwers, sig_throwers,
            mu_top=mu_top, log=True, full_output=True, llh_type=llh_type))
    # stuff.T_thresh = T_thresh
    # stuff.rate = pct
    # stuff.n_sigma = n_sigma
    save (stuff, 'sens.pickle')
            


                


