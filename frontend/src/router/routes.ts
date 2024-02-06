const drugDiscoveryRoutes = [
    {
        path: 'drug-discovery',
        name: 'DrugDiscovery',
        component: () => import('src/features/drug_discovery/ExperimentsView.vue')
    },
    {
        path: 'drug-discovery/experiment/:experimentId',
        component: () => import('src/features/drug_discovery/ExperimentNavigation.vue'),
        name: 'DrugDiscoveryExperiment',
        props: true,
        children: [
            {
                path: '/drug-discovery/experiment/:experimentId/upload-targets',
                component: () => import('../features/drug_discovery/components/UploadTargetsStep.vue'),
                name: 'UploadTargets',
                props: true
            },
            {
                path: '/drug-discovery/experiment/:experimentId/upload-ligands',
                component: () => import('../features/drug_discovery/components/UploadLigandsStep.vue'),
                name: 'UploadLigands',
                props: true
            },
            {
              path: '/drug-discovery/experiment/:experimentId/run-docking',
              component: () => import('../features/drug_discovery/components/RunDockingStep.vue'),
              name: 'RunDocking',
              props: true
            },
        ]
    }
];


const proteinDesignRoutes = [
    {
        path: 'protein-design',
        name: 'ProteinDesign',
        component: () => import('src/features/proteinDesign/ExperimentsView.vue'),
    },
    {
        path: 'protein-design/experiment/:experimentId',
        component: () => import('src/features/proteinDesign/ExperimentView.vue'),
        name: 'ProteinDesignExperiment',
        props: true
    },
];

const foldingRoutes = [
    {
        path: 'folding',
        name: 'Folding',
        component: () => import('src/features/aminoAcid/folding/ExperimentsView.vue'),
    },
    {
        path: 'folding/experiment/:experimentId',
        component: () => import('src/features/aminoAcid/folding/ExperimentView.vue'),
        name: 'FoldingExperiment',
        props: true
    },
];

const localisationRoutes = [
    {
        path: 'localisation',
        name: 'Localisation',
        component: () => import('src/features/aminoAcid/localisation/ExperimentsView.vue'),
    },
    {
        path: 'localisation/experiment/:experimentId',
        component: () => import('src/features/aminoAcid/localisation/ExperimentView.vue'),
        name: 'LocalisationExperiment',
        props: true
    },
];

const geneOntologyRoutes = [
    {
        path: 'gene-ontology',
        name: 'GeneOntology',
        component: () => import('src/features/aminoAcid/geneOntology/ExperimentsView.vue'),
    },
    {
        path: 'gene-ontology/experiment/:experimentId',
        component: () => import('src/features/aminoAcid/geneOntology/ExperimentView.vue'),
        name: 'GeneOntologyExperiment',
        props: true
    },
];

const solubilityRoutes = [
    {
        path: 'solubility',
        name: 'Solubility',
        component: () => import('src/features/aminoAcid/solubility/ExperimentsView.vue'),
    },
    {
        path: 'solubility/experiment/:experimentId',
        component: () => import('src/features/aminoAcid/solubility/ExperimentView.vue'),
        name: 'SolubilityExperiment',
        props: true
    },
];

const conformationsRoutes = [
    {
        path: 'conformations',
        name: 'Conformations',
        component: () => import("src/features/conformations/ExperimentsView.vue")
    },
    {
        path: 'conformations/experiment/:experimentId',
        component: () => import('src/features/conformations/ExperimentView.vue'),
        name: 'ConformationsExperiment',
        props: true
    }
]

const labsRoute = [
    {
        path: 'labs',
        name: 'Labs',
        component: () => import("src/components/LabsView.vue")
    },
];

const dt = [
    {
        path: 'drug-target2',
        name: 'DrugTarget2',
        component: () => import("src/features/drugTarget2/MainDrugTarget.vue")
    }
];

const routes = [
    {
        path: '/',
        name: 'Main',
        component: () => import('src/layouts/MainLayout.vue'),
        children: [
            ...dt,
            ...labsRoute,
            ...drugDiscoveryRoutes,
            ...proteinDesignRoutes,
            ...conformationsRoutes,
            ...localisationRoutes,
            ...solubilityRoutes,
            ...geneOntologyRoutes,
            ...foldingRoutes
        ]
    }
];

export default routes
