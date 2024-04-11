const drugDiscoveryRoutes = [
  {
    path: 'drug-discovery',
    name: 'Drug discovery experiments',
    component: () => import('src/features/drug_discovery/ExperimentsView.vue')
  },
  {
    path: 'drug-discovery/experiment/:experimentId',
    component: () => import('src/features/drug_discovery/ExperimentNavigation.vue'),
    name: 'Drug discovery',
    props: true,
    children: [
      {
        path: '/drug-discovery/experiment/:experimentId/upload-targets',
        component: () => import('../features/drug_discovery/components/UploadTargetsStep.vue'),
        name: 'Upload targets',
        props: true
      },
      {
        path: '/drug-discovery/experiment/:experimentId/upload-ligands',
        component: () => import('../features/drug_discovery/components/UploadLigandsStep.vue'),
        name: 'Upload ligands',
        props: true
      },
      {
        path: '/drug-discovery/experiment/:experimentId/run-docking',
        component: () => import('../features/drug_discovery/components/RunDockingStep.vue'),
        name: 'Run docking',
        props: true
      },
    ]
  }
];


const proteinDesignRoutes = [
  {
    path: 'protein-design',
    name: 'Protein design experiments',
    component: () => import('src/features/proteinDesign/ExperimentsView.vue'),
  },
  {
    path: 'protein-design/experiment/:experimentId',
    component: () => import('src/features/proteinDesign/ExperimentView.vue'),
    name: 'Protein design',
    props: true
  },
];

const foldingRoutes = [
  {
    path: 'folding',
    name: 'Folding experiments',
    component: () => import('src/features/aminoAcid/folding/ExperimentsView.vue'),
  },
  {
    path: 'folding/experiment/:experimentId',
    component: () => import('src/features/aminoAcid/folding/ExperimentView.vue'),
    name: 'Folding',
    props: true
  },
];

const localisationRoutes = [
  {
    path: 'localisation',
    name: 'Localisation experiments',
    component: () => import('src/features/aminoAcid/localisation/ExperimentsView.vue'),
  },
  {
    path: 'localisation/experiment/:experimentId',
    component: () => import('src/features/aminoAcid/localisation/ExperimentView.vue'),
    name: 'Localisation',
    props: true
  },
];

const geneOntologyRoutes = [
  {
    path: 'gene-ontology',
    name: 'Gene ontology experiments',
    component: () => import('src/features/aminoAcid/geneOntology/ExperimentsView.vue'),
  },
  {
    path: 'gene-ontology/experiment/:experimentId',
    component: () => import('src/features/aminoAcid/geneOntology/ExperimentView.vue'),
    name: 'Gene ontology',
    props: true
  },
];

const solubilityRoutes = [
  {
    path: 'solubility',
    name: 'Solubility experiments',
    component: () => import('src/features/aminoAcid/solubility/ExperimentsView.vue'),
  },
  {
    path: 'solubility/experiment/:experimentId',
    component: () => import('src/features/aminoAcid/solubility/ExperimentView.vue'),
    name: 'Solubility',
    props: true
  },
];

const conformationsRoutes = [
  {
    path: 'conformations',
    name: 'Conformations experiments',
    component: () => import("src/features/conformations/ExperimentsView.vue")
  },
  {
    path: 'conformations/experiment/:experimentId',
    component: () => import('src/features/conformations/ExperimentView.vue'),
    name: 'Conformations',
    props: true
  }
]

const smallMoleculeDesignRoutes = [
  {
    path: 'small-molecules',
    name: 'Small molecules design experiments',
    component: () => import("src/features/smallMoleculeDesign/ExperimentsView.vue")
  },
  {
    path: 'small-molecules/experiment/:experimentId',
    component: () => import('src/features/smallMoleculeDesign/ExperimentView.vue'),
    name: 'Small molecules design',
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

const routes = [
  {
    path: '/',
    name: 'Main',
    component: () => import('src/layouts/MainLayout.vue'),
    children: [
      ...labsRoute,
      ...drugDiscoveryRoutes,
      ...proteinDesignRoutes,
      ...conformationsRoutes,
      ...localisationRoutes,
      ...solubilityRoutes,
      ...geneOntologyRoutes,
      ...foldingRoutes,
      ...smallMoleculeDesignRoutes
    ]
  }
];

export default routes
