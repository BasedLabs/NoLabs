<template>
  <q-page class="bg-black q-pl-md">
    <BioBuddyChat v-if="bioBuddyEnabled" :experiment-id="experiment.experimentId" />
    <div class="row no-wrap items-center">
      <button @click="generateWorkflow">Generate Workflow</button>
      <q-btn-dropdown color="primary" label="Add Node" icon="add" dense persistent>
        <q-list>
          <q-item v-for="option in componentOptions" :key="option.value" clickable v-close-popup
            @click="addComponent(option)">
            <q-item-section>
              <q-item-label>{{ option.label }}</q-item-label>
              <q-item-label caption>{{ option.description }}</q-item-label>
            </q-item-section>
          </q-item>
        </q-list>
      </q-btn-dropdown>
      <q-space />
      <q-btn class="q-ma-sm" color="green" icon="not_started">Run workflow</q-btn>
    </div>

    <div class="map-container">
      <VueFlow class="workflow" :elevate-edges-on-select="true" v-if="elements" :nodes="elements.nodes" @nodeDragStop="onNodeDragStopHandler"
        :edges="elements.edges" @connect="onConnect" fit-view-on-init>
        <template #node-protein-list="{ id }">
          <ProteinListNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-ligand-list="{ id }">
          <LigandListNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-esmfold="{ id }">
          <EsmFoldNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-diffdock="{ id }">
          <DiffDockNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
        </template>
        <template #node-rfdiffusion="{ id }">
          <RfDiffusionNode :nodeId="id" :inputs="elements.nodes.find(n => n.id === id)?.data.inputs"
            :outputs="elements.nodes.find(n => n.id === id)?.data.outputs" :onDeleteNode="onDeleteNode"
            :onOpenSettings="openSettings" :onOpenDialog="openNodeDialog" />
        </template>
      </VueFlow>
    </div>
  </q-page>

  <q-drawer v-model="sideMenuOpen" v-show="sideMenuOpen" bordered content-style="background-color: white;" side="right">
    <!-- Close button -->
    <q-btn @click="closeSideMenu" class="q-ma-md" flat round dense icon="close" />

    <!-- Side menu content -->
    <q-card>
      <q-card-section>
        <q-item>
          <q-item-section>
            <q-item-label v-if="selectedNode">{{ selectedNode.label }}</q-item-label>
            <q-item-label v-if="selectedNode" caption>{{ selectedNode.description }}</q-item-label>
          </q-item-section>
        </q-item>
      </q-card-section>
    </q-card>
  </q-drawer>

  <q-dialog v-model="modalOpen" persistent>
    <q-card style="min-width: 70vw; min-height: 70vh;">
      <q-card-actions align="right">
        <q-btn flat round dense icon="close" v-close-popup @click="closeModal" />
      </q-card-actions>
      <q-card-section>
        <ProteinListNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'protein-list'"
          :experiment-id="experiment.experimentId" />
        <LigandListNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'ligand-list'"
          :experiment-id="experiment.experimentId" />
        <EsmFoldNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'esmfold'"
          :experiment-id="experiment.experimentId" />
        <DiffDockNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'diffdock'"
          :experiment-id="experiment.experimentId" />
        <RfDiffusionNodeContent v-if="experiment.experimentId && selectedNode && selectedNode.type === 'rfdiffusion'"
          :experiment-id="experiment.experimentId" />
      </q-card-section>
    </q-card>
  </q-dialog>
</template>

<script lang="ts">
import '@vue-flow/core/dist/style.css';
import '@vue-flow/core/dist/theme-default.css';

import BioBuddyChat from "src/features/biobuddy/BioBuddyChat.vue";
import ProteinListNode from "./components/workflow/ProteinListNode.vue";
import ProteinListNodeContent from "./components/workflow/ProteinListNodeContent.vue";
import LigandListNode from "./components/workflow/LigandListNode.vue";
import LigandListNodeContent from "./components/workflow/LigandListNodeContent.vue";
import { useDrugDiscoveryStore } from "./storage";
import { defineComponent } from "vue";
import { checkBioBuddyEnabled } from "../biobuddy/api";
import { Edge, Position, Node as FlowNode } from '@vue-flow/core';
import { VueFlow } from '@vue-flow/core';
import EsmFoldNode from "./components/workflow/EsmFoldNode.vue";
import EsmFoldNodeContent from "./components/workflow/EsmFoldNodeContent.vue";
import DiffDockNode from "./components/workflow/DiffDockNode.vue";
import DiffDockNodeContent from "./components/workflow/DiffDockNodeContent.vue";
import RfDiffusionNode from "./components/workflow/RfDiffusionNode.vue";
import RfDiffusionNodeContent from "./components/workflow/RfDiffusionNodeContent.vue";

// Define custom Node type
interface Node extends FlowNode {
  id: string;
  name: string;
  type: string;
  inputs: string[];
  outputs: string[];
  description: string;
}

export default defineComponent({
  name: "ExperimentNavigation",
  components: {
    BioBuddyChat,
    VueFlow,
    ProteinListNode,
    ProteinListNodeContent,
    LigandListNode,
    LigandListNodeContent,
    EsmFoldNode,
    EsmFoldNodeContent,
    DiffDockNode,
    DiffDockNodeContent,
    RfDiffusionNode,
    RfDiffusionNodeContent
  },
  data() {
    return {
      step: 1, // This can be a reactive property based on the current route if needed
      experiment: {
        experimentId: null as string | null,
        metadata: null
      },
      sideMenuOpen: false,
      modalOpen: false,
      selectedNode: null as Node | null,
      selectedComponent: null as string | null,
      specialNodeProps: [],
      splitterModel: 20,
      bioBuddyEnabled: false,
      componentOptions: [
        { label: 'Protein List Node', value: 'protein-list', description: 'Downloads target protein data from RCSB PDB in FASTA format.' },
        { label: 'Ligand List Node', value: 'ligand-list', description: 'Downloads ligands from ChEMBL in SMILES format.' },
        { label: 'Input Node', value: 'input', description: 'Represents an input node.' },
        { label: 'Output Node', value: 'output', description: 'Represents an output node.' },
      ],
      elements: {
        nodes: [] as Node[],
        edges: [] as Edge[]
      },
      jsonData: {
        "workflow": {
          "name": "Protein-Ligand Interaction and Protein Generation Workflow",
          "nodes": [
            {
              "id": "1",
              "name": "DownloadFromRCSB",
              "type": "component",
              "inputs": [],
              "outputs": ["pdb_file"],
              "description": "Downloads target protein data from RCSB PDB."
            },
            {
              "id": "2",
              "name": "DownloadFromChembl",
              "type": "component",
              "inputs": [],
              "outputs": ["smiles_string"],
              "description": "Downloads ligands from ChEMBL in SMILES format."
            },
            {
              "id": "3",
              "name": "DockingProteinOnLigand",
              "type": "diffdock",
              "inputs": ["pdb_file", "smiles_string"],
              "outputs": ["pdb_file"],
              "description": "Performs docking of protein on ligand."
            },
            {
              "id": "4",
              "name": "GenerateProteins",
              "type": "rfdiffusion",
              "inputs": [],
              "outputs": ["generated_fasta_file"],
              "description": "Generates new protein sequences."
            },
            {
              "id": "5",
              "name": "RunFolding",
              "type": "esmfold",
              "inputs": ["generated_fasta_file"],
              "outputs": ["generated_pdb_file"],
              "description": "Folds generated protein sequences into 3D structures."
            },
            {
              "id": "6",
              "name": "CalculateBlast",
              "type": "component",
              "inputs": ["protein_complex_pdb_file"],
              "outputs": ["similarity_scores"],
              "description": "Calculates similarity scores between protein complexes."
            }
          ],
          "edges": [
            { "from": "1-output-pdb_file", "to": "3-input-pdb_file" },
            { "from": "2-output-smiles_string", "to": "3-input-smiles_string" },
            { "from": "4-output-generated_fasta_file", "to": "5-input-generated_fasta_file" },
            { "from": "5-output-generated_pdb_file", "to": "3-input-generated_pdb_file" },
            { "from": "1-output-pdb_file", "to": "6-input-pdb_file" }
          ]
        }
      }
    };
  },
  async mounted() {
    const store = useDrugDiscoveryStore();
    this.setStepBasedOnRoute();
    this.experiment.experimentId = this.$route.params.experimentId as string;
    this.experiment.metadata = await store.getExperimentMetaData(this.experiment.experimentId);
    try {
      const response = await checkBioBuddyEnabled();
      this.bioBuddyEnabled = response.enabled;
    } catch (error) {
      console.error('Error checking BioBuddy enabled status:', error);
      this.bioBuddyEnabled = false;
    }
  },
  methods: {
    setStepBasedOnRoute() {
      switch (this.$route.name) {
        case 'Upload targets':
          this.step = 1;
          break;
        case 'Upload ligands':
          this.step = 2;
          break;
        case 'Run docking':
          this.step = 3;
          break;
        default:
          this.step = 1; // Default step if the route name doesn't match
      }
    },
    generateWorkflow() {
    const nodes: Node[] = [];
    const edges: Edge[] = [];

    // Helper function to get the depth of a node in the graph
    const getNodeDepth = (nodeId: string): number => {
      const incomingEdges = this.jsonData.workflow.edges.filter(edge => edge.to.split('-')[0] === nodeId);
      if (incomingEdges.length === 0) return 0; // Node has no incoming edges, depth is 0
      const parentDepths = incomingEdges.map(edge => getNodeDepth(edge.from.split('-')[0]));
      return Math.max(...parentDepths) + 1; // Depth is the maximum parent depth plus 1
    };

    // Calculate depths for all nodes
    const depths: { [key: string]: number } = {};
    this.jsonData.workflow.nodes.forEach(node => {
      depths[node.id] = getNodeDepth(node.id);
    });

    // Group nodes by depth
    const nodesByDepth: { [key: number]: Node[] } = {};
    this.jsonData.workflow.nodes.forEach(node => {
      const depth = depths[node.id];
      if (!nodesByDepth[depth]) nodesByDepth[depth] = [];
      nodesByDepth[depth].push(node);
    });

    // Calculate positions for each node
    let currentX = 100;
    const layerWidth = 500; // Width difference between layers
    const layerHeight = 300; // Height difference between nodes within a layer
    Object.values(nodesByDepth).forEach((nodesAtDepth: Node[], depth: number) => {
      const numNodes = nodesAtDepth.length;
      const startY = 100 + (numNodes * layerHeight) / 2;

      nodesAtDepth.forEach((node, index) => {
        const position = { x: currentX, y: startY - index * layerHeight };
        const nodeData: Node = {
          id: node.id,
          name: node.name,
          type: node.type,
          data: {
            description: node.description,
            inputs: node.inputs,
            outputs: node.outputs,
            draggable: false
          },
          position
        };

        nodes.push(nodeData);
      });

      // Update currentX for the next layer
      currentX += layerWidth;
    });

    // Add edges for nested nodes
    this.jsonData.workflow.edges.forEach(edge => {
      edges.push({
        id: `e${edge.from}-to-${edge.to}`,
        source: edge.from.split('-')[0],
        target: edge.to.split('-')[0],
        sourceHandle: edge.from,
        targetHandle: edge.to,
        style: { zIndex: 1 } // Ensure edges are on top
      });
    });

    this.elements = { nodes, edges };
    },
    onConnect(params: { source: string, target: string }) {
      const newEdge: Edge = {
        id: `e${params.source}-${params.target}`,
        source: params.source,
        target: params.target,
      };
      this.elements.edges.push(newEdge);
    },
    onDeleteNode(nodeId: string) {
      const index = this.elements.nodes.findIndex(node => node.id === nodeId);
      if (index !== -1) {
        this.elements.nodes.splice(index, 1);
      }
    },
    openSettings(nodeId: string) {
      const node = this.elements.nodes.find(node => node.id === nodeId);
      if (node) {
        this.selectedNode = node;
        this.sideMenuOpen = true;
      }
    },
    closeSideMenu() {
      this.sideMenuOpen = false;
    },
    closeModal() {
      this.modalOpen = false;
    },
    onNodeDragStopHandler(event: { node: Node }) {
      const id = event.node.id;
      const newPosition = event.node.position;

      // Find the node with the given id and update its position
      const nodeToUpdate = this.elements.nodes.find((node: Node) => node.id === id);
      if (nodeToUpdate) {
        nodeToUpdate.position = newPosition;
      }
    },
    async onExperimentNameChange(newExperimentName: string) {
      const store = useDrugDiscoveryStore();
      await store.changeExperimentName(this.experiment.experimentId as string, newExperimentName);
    },
    openNodeDialog(nodeId: string) {
      const node = this.elements.nodes.find(node => node.id === nodeId);
      if (node) {
        this.selectedNode = node;
        this.modalOpen = true;
      }
    }
  }
})
</script>

<style scoped>
body {
  overflow-x: hidden;
}

.map-container {
  width: 100vw;
  height: 100vh;
}

.workflow {
  height: 100%;
  width: 100%;
}

.vue-flow__edge {
  z-index: 2;
  position: relative;
}

.vue-flow__node {
  z-index: 1;
  position: relative;
}
</style>
