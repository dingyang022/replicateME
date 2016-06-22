from numpy import array
from scipy.sparse import dok_matrix
import cobra
from cobra import Model, DictList
from cobra import Metabolite as Component

from reactionfiles import MEReaction,StoichiometricData
from sympy import Symbol
import pandas
from six import iteritems


mu = Symbol("mu", positive=True)

class MEComponent(Component):

    def __init__(self, id):
        Component.__init__(self, id)
        pass

    def remove_from_MEmodel(self, method='subtractive'):
        try:
            self._model.process_data.remove(self.id)
            self._model.complex_data.remove(self.id)
        except:
            pass
        if method == 'subtractive':
		self.remove_from_model(method=method)

class TranscribedGene(MEComponent):

    def __init__(self, id):
        MEComponent.__init__(self, id)
        self.left_pos = None
        self.right_pos = None
        self.strand = None
        self.RNA_type = ''
        self.nucleotide_sequence = ''

    @property
    def nucleotide_count(self):
        seq = self.nucleotide_sequence
        counts = {i: seq.count(i) for i in ("A", "T", "G", "C")}
        monophosphate_counts = {dogma.transcription_table[k].replace("tp_c",
                                                                     "mp_c"): v
                                for k, v in iteritems(counts)}
        return monophosphate_counts

    @property
    def mass(self):
	return compute_RNA_mass(self.nucleotide_sequence)

class Complex(MEComponent):
    @property
    def metabolic_reactions(self):
        """read-only link to MetabolicReactions"""
        reaction_list = []
        for reaction in self.reactions:
            if reaction.__class__.__name__ == 'MetabolicReaction':
                reaction_list.append(reaction)
	return reaction_list

def add_m_model_content(me_model, m_model, complex_metabolite_ids=[]):
    """
    Add metabolite and reaction attributes to me_model from m_model. Also
    creates StoichiometricData objects for each reaction in m_model, and adds
    reactions directly to me_model if they are exchanges or demands.
    Args:
        me_model: cobra.model.MEModel
            The MEModel object to which the content will be added
        m_model: cobra.model
            The m_model which will act as the source of metabolic content for
            MEModel
        complex_metabolite_ids: list
            List of complexes which are 'metabolites' in the m-model reaction
            matrix, but should be treated as complexes
    """
    for met in m_model.metabolites:
        if met.id in complex_metabolite_ids:
            new_met = Complex(met.id)
        elif met.id.startswith("RNA"):
            new_met = TranscribedGene(met.id)
        else:
            new_met = MEComponent(met.id)
        new_met.name = met.name
        new_met.formula = met.formula
        new_met.compartment = met.compartment
        new_met.charge = met.charge
        new_met.annotation = met.annotation
        new_met.notes = met.notes
        me_model.add_metabolites(new_met)

    for reaction in m_model.reactions:
        if reaction.id.startswith("EX_") or reaction.id.startswith("DM_"):
            new_reaction = cobra.Reaction(reaction.id)
            me_model.add_reaction(new_reaction)
            new_reaction.lower_bound = reaction.lower_bound
            new_reaction.upper_bound = reaction.upper_bound
            for met, stoichiometry in iteritems(reaction.metabolites):
                new_reaction.add_metabolites(
                    {me_model.metabolites.get_by_id(met.id): stoichiometry})

        else:
            reaction_data = StoichiometricData(reaction.id, me_model)
            reaction_data.lower_bound = reaction.lower_bound
            reaction_data.upper_bound = reaction.upper_bound
            reaction_data._stoichiometry = {k.id: v for k, v in iteritems(reaction.metabolites)}


def get_reaction_matrix_dict():
    reaction_matrix = open('./raw_data/reaction_matrix.txt', 'r')
    # These metabolites are mistakenly labeled as NoCompartment when they
    # should really be in the cytosol.
    move_to_cytosol = {'adp', 'atp', 'h', 'pi', '2tpr3dpcoa', 'dpm', 'fe2',
                       'dad__5', 'met__L', 'tl'}
    ME_reaction_dict = {}
    for line in reaction_matrix:
        line = line.strip()
        if line.startswith("#") or len(line) == 0:
            continue
        rxn, met, comp, count = line.split('\t')
        rxn = rxn.replace('DASH', '')
        met = met.replace('DASH', '')
        # use compartment to append appropriate suffix
        if comp == 'Cytosol':
            met += '_c'
        elif comp == 'Periplasm':
            met += '_p'
        elif comp == 'Extra-organism':
            met += '_e'
        # some mistakenly annotated as no compartment
        elif comp == 'No_Compartment' and met in move_to_cytosol:
            met += '_c'
        if rxn not in ME_reaction_dict:
            ME_reaction_dict[rxn] = {}
        ME_reaction_dict[rxn][met] = float(count)
    reaction_matrix.close()
    for rxn_id in ["PFL_act", "hemeD_synthesis", "23bpg_generation"]:
        ME_reaction_dict[rxn_id] = \
            {k + "_c": v for k, v in iteritems(ME_reaction_dict[rxn_id])}
    return ME_reaction_dict


def get_reaction_info_frame():
    return pandas.read_csv("./raw_data/reactions.txt",
                           delimiter="\t", index_col=0)

def fix_id(id_str):
    return id_str.replace("_DASH_", "__")

def get_m_model():
    m = cobra.Model("e_coli_ME_M_portion")
    m.compartments = {"p": "Periplasm", "e": "Extra-organism", "c": "Cytosol"}
    compartment_lookup = {v: k for k, v in m.compartments.items()}

    met_info = pandas.read_csv("./raw_data/metabolites.txt",
                               delimiter="\t", header=None, index_col=0,
                               names=["id", "name", "formula", "compartment",
                                      "data_source"])

    for met_id in met_info.index:
        fixed_id = fix_id(met_id)
        for compartment in met_info.compartment[met_id].split("AND"):
            compartment = compartment.strip()
            if compartment == "No_Compartment":
                print("Assigned %s to c" % met_id)
                compartment = m.compartments["c"]
            new_met = cobra.Metabolite(
                fixed_id + "_" + compartment_lookup[compartment])
            new_met.name = met_info.name[met_id]
            new_met.formula = met_info.formula[met_id]
            m.add_metabolites(new_met)

    rxn_info = get_reaction_info_frame()
    rxn_dict = get_reaction_matrix_dict()
    for rxn_id in rxn_info.index:
        reaction = cobra.Reaction(rxn_id)
        reaction.name = rxn_info.description[rxn_id]
        for met_id, amount in rxn_dict[rxn_id].items():
            try:
                metabolite = m.metabolites.get_by_id(met_id)
            except KeyError:
                metabolite = cobra.Metabolite(met_id)
            reaction.add_metabolites({metabolite: amount})
        reaction.lower_bound = \
            -1000. if rxn_info.is_reversible[rxn_id] else 0.
        reaction.upper_bound = 1000.
        if rxn_info.is_spontaneous[rxn_id]:
            reaction.gene_reaction_rule = "s0001"
        m.add_reaction(reaction)

    sources_sinks = pandas.read_csv(
        "./raw_data/reaction_matrix_sources_and_sinks.txt",
        delimiter="\t", header=None, names=["rxn_id", "met_id", "compartment",
                                            "stoic"], index_col=1)

    source_amounts = pandas.read_csv("./raw_data/exchange_bounds.txt",
                                     delimiter="\t", index_col=0,
                                     names=["met_id", "amount"])

    sources_sinks.index = [fix_id(i) for i in sources_sinks.index]
    source_amounts.index = [fix_id(i) for i in source_amounts.index]

    for met in sources_sinks.index:
        met_id = met + "_" + compartment_lookup[sources_sinks.compartment[met]]
        # EX_ or DM_ + met_id
        reaction_id = sources_sinks.rxn_id[met][:3] + met_id
        reaction = cobra.Reaction(reaction_id)
        m.add_reaction(reaction)
        reaction.add_metabolites({m.metabolites.get_by_id(met_id): -1})
        # set bounds on exchanges
        if reaction.id.startswith("EX_") and met in source_amounts.index:
            reaction.lower_bound = -source_amounts.amount[met]
    return m

class MEmodel(Model):
    def __init__(self, *args):
        Model.__init__(self, *args)
        self.global_info = {}
        self.stoichiometric_data = DictList()
        self.complex_data = DictList()
        self.modification_data = DictList()
        self.translation_data = DictList()
        self.transcription_data = DictList()
        self.generic_data = DictList()
        self.tRNA_data = DictList()
        self.translocation_data = DictList()
        self.posttranslation_data = DictList()
        self.subreaction_data = DictList()
        self.process_data = DictList()
        # create the biomass/dilution constraint
        self._biomass = Constraint("biomass")
        self._biomass_dilution = SummaryVariable("biomass_dilution")
        self._biomass_dilution.add_metabolites({self._biomass: -1})
        self.add_reaction(self._biomass_dilution)
        self._biomass_dilution.upper_bound = mu
        self._biomass_dilution.lower_bound = mu
        # Unmodeled protein is handled by converting protein_biomass to
        # biomass, and requiring production of the appropriate amount of dummy
        # protein
        self._unmodeled_protein_fraction = None
        self._protein_biomass = Constraint("protein_biomass")
        self._protein_biomass_dilution = SummaryVariable("protein_biomass_dilution")
        self._protein_biomass_dilution.add_metabolites({
            self._protein_biomass: -1,
            self._biomass: 1,
        })
        self._mRNA_biomass = Constraint("mRNA_biomass")
        self._mRNA_biomass_dilution = SummaryVariable("mRNA_biomass_dilution")
        self._mRNA_biomass_dilution.add_metabolites({
            self._mRNA_biomass: -1,
            self._biomass: 1,
        })
        self._tRNA_biomass = Constraint("tRNA_biomass")
        self._tRNA_biomass_dilution = SummaryVariable("tRNA_biomass_dilution")
        self._tRNA_biomass_dilution.add_metabolites({
            self._tRNA_biomass: -1,
            self._biomass: 1,
        })
        self._rRNA_biomass = Constraint("rRNA_biomass")
        self._rRNA_biomass_dilution = SummaryVariable("rRNA_biomass_dilution")
        self._rRNA_biomass_dilution.add_metabolites({
            self._rRNA_biomass: -1,
            self._biomass: 1,
        })

        self._ncRNA_biomass = Constraint("ncRNA_biomass")
        self._ncRNA_biomass_dilution = SummaryVariable("ncRNA_biomass_dilution")
        self._ncRNA_biomass_dilution.add_metabolites({
            self._ncRNA_biomass: -1,
            self._biomass: 1,
        })
        self.add_reactions((self._protein_biomass_dilution,
                            self._mRNA_biomass_dilution,
                            self._tRNA_biomass_dilution,
                            self._rRNA_biomass_dilution,
                            self._ncRNA_biomass_dilution))

        self._DNA_biomass = Constraint("DNA_biomass")
        self._DNA_biomass_dilution = SummaryVariable("DNA_biomass_dilution")
        self._DNA_biomass_dilution.add_metabolites({
            self._DNA_biomass: -1e-3,
            self._biomass: 1e-3,
        })
        self._DNA_biomass_dilution.lower_bound = mu
        self._DNA_biomass_dilution.upper_bound = mu

    @property
    def unmodeled_protein(self):
        return self.metabolites.get_by_id("protein_dummy")

    @property
    def unmodeled_protein_fraction(self):
        return self._unmodeled_protein_fraction

    @unmodeled_protein_fraction.setter
    def unmodeled_protein_fraction(self, value):
        # proportion = value / (1 - value)
        # see the Biomass_formulations for an explanation
        amount = value / self.unmodeled_protein.mass
        self._protein_biomass_dilution.add_metabolites(
                {self.unmodeled_protein: -amount}, combine=False)
        self._protein_biomass_dilution.add_metabolites(
            {self._biomass: 1+value}, combine=False)
        self._unmodeled_protein_fraction = value

    def get_metabolic_flux(self, solution=None):
        """extract the flux state for metabolic reactions"""
        if solution is None:
            solution = self.solution
        if solution.status != "optimal":
            raise ValueError("solution status '%s' is not 'optimal'" %
                             solution.status)
        flux_dict = {r.id: 0 for r in self.stoichiometric_data}
        for reaction in self.reactions:
            if isinstance(reaction, MetabolicReaction):
                m_reaction_id = reaction.stoichiometric_data.id
                if reaction.reverse:
                    flux_dict[m_reaction_id] -= solution.x_dict[reaction.id]
                else:
                    flux_dict[m_reaction_id] += solution.x_dict[reaction.id]
            elif reaction.id.startswith("EX_") or reaction.id.startswith("DM"):
                flux_dict[reaction.id] = solution.x_dict[reaction.id]
        return flux_dict

    def get_transcription_flux(self, solution=None):
        """extract the transcription flux state"""
        if solution is None:
            solution = self.solution
        if solution.status != "optimal":
            raise ValueError("solution status '%s' is not 'optimal'" %
                             solution.status)
        flux_dict = {}
        for reaction in self.reactions:
            if isinstance(reaction, TranscriptionReaction):
                for rna_id in reaction.transcription_data.RNA_products:
                    locus_id = rna_id.replace("RNA_", "", 1)
                    if locus_id not in flux_dict:
                        flux_dict[locus_id] = 0
                    flux_dict[locus_id] += solution.x_dict[reaction.id]
        return flux_dict

    def get_translation_flux(self, solution=None):
        """extract the translation flux state"""
        if solution is None:
            solution = self.solution
        if solution.status != "optimal":
            raise ValueError("solution status '%s' is not 'optimal'" %
                             solution.status)
        flux_dict = {r.id: 0 for r in self.translation_data}
        for reaction in self.reactions:
            if isinstance(reaction, TranslationReaction):
                protein_id = reaction.translation_data.id
                flux_dict[protein_id] += solution.x_dict[reaction.id]
        return flux_dict

    def construct_S(self, growth_rate):
        """build the stoichiometric matrix at a specific growth rate"""
        # intialize to 0
        S = dok_matrix((len(self.metabolites), len(self.reactions)))
        # populate with stoichiometry
        for i, r in enumerate(self.reactions):
            for met, value in r._metabolites.iteritems():
                met_index = self.metabolites.index(met)
                if hasattr(value, "subs"):
                    S[met_index, i] = float(value.subs(mu, growth_rate))
                else:
                    S[met_index, i] = float(value)
        return S

    def construct_attribute_vector(self, attr_name, growth_rate):
        """build a vector of a reaction attribute at a specific growth rate

        Mainly used for upper and lower bounds"""
        return array([float(value.subs(mu, growth_rate))
                      if hasattr(value, "subs") else float(value)
                      for value in self.reactions.list_attr(attr_name)])

    def compute_solution_error(self, solution=None):
        errors = {}
        if solution is None:
            solution = self.solution
        S = self.construct_S(solution.f)
        lb = self.construct_attribute_vector("lower_bound", solution.f)
        ub = self.construct_attribute_vector("upper_bound", solution.f)
        x = array(solution.x)
        err = abs(S * x)
        errors["max_error"] = err.max()
        errors["sum_error"] = err.sum()
        ub_err = min(ub - x)
        errors["upper_bound_error"] = abs(ub_err) if ub_err < 0 else 0
        lb_err = min(x - lb)
        errors["lower_bound_error"] = abs(lb_err) if lb_err < 0 else 0
        return errors

    def update(self):
        """updates all component reactions"""
        for r in self.reactions:
            if hasattr(r, "update"):
                r.update()

    def prune(self):
        """remove all unused metabolites and reactions

        This should be run after the model is fully built. It will be
        difficult to add new content to the model once this has been run.

        """
        complex_data_list = [i.id for i in self.complex_data]
        for c_d in complex_data_list:
            c = self.complex_data.get_by_id(c_d)
            cplx = c.complex
            if len(cplx.reactions) == 1:
                list(cplx.reactions)[0].delete(remove_orphans=True)
                self.complex_data.remove(self.complex_data.get_by_id(c_d))

        for p in self.metabolites.query(re.compile('^protein_')):
            if isinstance(p, ProcessedProtein):
                delete = True
                for rxn in p._reaction:
                    try:
                        if p in rxn.reactants:
                            delete = False
                    except Exception as e:
                        print(rxn)
                        raise e
                if delete:
                    while len(p._reaction) > 0:
                        list(p._reaction)[0].delete(remove_orphans=True)
                        for data in self.posttranslation_data.query(p.id):
                            self.posttranslation_data.remove(data.id)

        for p in self.metabolites.query(re.compile('^protein_')):
            if isinstance(p, TranslatedGene):
                delete = True
                for rxn in p._reaction:
                    try:
                        if p in rxn.reactants and not rxn.id.startswith('degradation'):
                            delete = False
                    except Exception as e:
                        print(rxn)
                        raise e
                if delete:
                    while len(p._reaction) > 0:
                        list(p._reaction)[0].delete(remove_orphans=True)
                        p_id = p.id.replace('protein_', '')
                        for data in self.translation_data.query(p_id):
                            self.translation_data.remove(data.id)


        removed_RNA = set()
        for m in list(self.metabolites.query(re.compile("^RNA_"))):
            delete = True
            for rxn in m._reaction:
                if m in rxn.reactants and not rxn.id.startswith('DM_'):
                    delete = False
            if delete:
                try:
                    self.reactions.get_by_id('DM_' + m.id).remove_from_model(
                            remove_orphans=True)
                    if m in self.metabolites:
                        m.remove_from_model(method='subtractive')
                except KeyError:
                    pass
                else:
                    removed_RNA.add(m.id)

        for t in self.reactions.query('transcription_TU'):
            delete = True
            for product in t.products:
                if isinstance(product, TranscribedGene):
                    delete = False
            t_process_id = t.id.replace('transcription_', '')
            if delete:
                t.remove_from_model(remove_orphans=True)
                self.transcription_data.remove(t_process_id)
            else:
                # gets rid of the removed RNA from the products
                self.transcription_data.get_by_id(
                    t_process_id).RNA_products.difference_update(removed_RNA)

    def remove_genes_from_model(self, gene_list):
        for gene in gene_list:
            self.metabolites.get_by_id('RNA_'+gene).remove_from_model(method='subtractive')
            protein = self.metabolites.get_by_id('protein_'+gene)
            for cplx in protein.complexes:
                print cplx
                for rxn in cplx.metabolic_reactions:
                    try:
                        self.stoichiometric_data.remove(rxn.id.split('_')[0])
                    except ValueError:
                        pass
                    rxn.remove_from_model()

            protein.remove_from_model(method='destructive')

        # Remove all transcription reactions that now do not form a used
        # transcript
        for t in self.reactions.query('transcription_TU'):
            delete = True
            for product in t.products:
                if isinstance(product, TranscribedGene):
                    delete = False
            if delete:
                t.remove_from_model(remove_orphans=True)
                t_process_id = t.id.replace('transcription_', '')
                self.transcription_data.remove(t_process_id)

    def get_biomass_composition(self, solution=None):
        if solution is None:
            solution = self.solution
        biomass_composition = defaultdict(float)
        for met, stoich in self._protein_biomass_dilution.metabolites.items():
            if abs(stoich) < 1:
                weight = self.unmodeled_protein.mass
                biomass_composition['Unmodeled Protein'] = \
                    solution.x_dict['protein_biomass_dilution'] * \
                    abs(stoich) * weight
        biomass_composition['Protein'] = \
            solution.x_dict['protein_biomass_dilution']
        biomass_composition['tRNA'] = \
            solution.x_dict['tRNA_biomass_dilution']
        biomass_composition['mRNA'] = \
            solution.x_dict['mRNA_biomass_dilution']
        biomass_composition['ncRNA'] = \
            solution.x_dict['ncRNA_biomass_dilution']
        biomass_composition['rRNA'] = \
            solution.x_dict['rRNA_biomass_dilution']
        biomass_composition['Other'] = \
            solution.x_dict['biomass_component_dilution']

        return biomass_composition

    def RNA_to_protein_ratio(self, solution=None):
        composition = self.get_biomass_composition(solution=solution)
        RNA_to_protein = (composition['mRNA'] + composition['tRNA'] +
                          composition['rRNA'] + composition['ncRNA']) / \
                         (composition['Protein'] +
                          composition['Unmodeled Protein'])
        return RNA_to_protein

    def get_RNA_fractions_dict(self, solution=None):
        RNA_fractions = {}
        composition = self.get_biomass_composition(solution=solution)

        tRNA_to_RNA = (composition['tRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['tRNA'] = tRNA_to_RNA

        rRNA_to_RNA = (composition['rRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['rRNA'] = rRNA_to_RNA

        mRNA_to_RNA = (composition['mRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['mRNA'] = mRNA_to_RNA

        ncRNA_to_RNA = (composition['ncRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
        RNA_fractions['ncRNA'] = ncRNA_to_RNA

        return RNA_fractions

    def make_biomass_composition_piechart(self, solution=None):
        try:
            import brewer2mpl
        except ImportError:
            color_map = None
        else:
            color_map = brewer2mpl.wesanderson.Zissou.mpl_colormap

        try:
            import pandas
        except ImportError:
            raise Exception("Pandas must be installed to get biomass piechart")

        if solution is None:
            solution = self.solution

        summary = {}
        summary['Biomass composition'] = \
            self.get_biomass_composition(solution=solution)
        frame = pandas.DataFrame.from_dict(summary) / solution.f


        print 'Total biomass sum =', frame.sum().values[0]
        return frame.plot(kind='pie', subplots=True, legend=None, colormap=color_map)

