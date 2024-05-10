<template>
  <q-page class="bg-black q-pl-md q-pr-md">
    <BioBuddyChat v-if="bioBuddyEnabled" :experiment-id=" this.$route.params.experimentId" />
    <button @click="generateWorkflow">Generate Workflow</button>
    <q-btn-dropdown color="primary" label="Add Node" icon="add" dense persistent>
      <q-list>
        <q-item v-for="option in componentOptions" :key="option.value" clickable v-close-popup @click="addComponent(option)">
          <q-item-section>
            <q-item-label>{{ option.label }}</q-item-label>
            <q-item-label caption>{{ option.description }}</q-item-label>
          </q-item-section>
        </q-item>
      </q-list>
    </q-btn-dropdown>
    <div class="map-container">
      <VueFlow class="workflow"
               v-if="elements"
               :nodes="elements.nodes"
               @nodeDragStop="onNodeDragStopHandler"
               :edges="elements.edges"
               @connect="onConnect"
               fit-view-on-init>
        <template #node-protein-list>
          <ProteinListNode />
        </template>
        <template #node-ligand-list="{ id }">
          <LigandListNode :nodeId="id" :onDeleteNode="onDeleteNode" />
        </template>
      </VueFlow>
    </div>
  </q-page>

</template>

<script lang="ts">
import '@vue-flow/core/dist/style.css';
import '@vue-flow/core/dist/theme-default.css';

import BioBuddyChat from "src/features/biobuddy/BioBuddyChat.vue";
import {useDrugDiscoveryStore} from "./storage";
import {defineComponent} from "vue";
import {checkBioBuddyEnabled} from "../biobuddy/api";
import { ref } from 'vue';
import {Edge, Elements, Position, updateEdge, useVueFlow} from '@vue-flow/core';
import { VueFlow } from '@vue-flow/core';
import ProteinListNode from "./components/workflow/ProteinListNode.vue";
import LigandListNode from "./components/workflow/LigandListNode.vue";

export default defineComponent({
  name: "ExperimentNavigation",
  components: {BioBuddyChat, VueFlow, ProteinListNode, LigandListNode},
  data() {
    return {
      step: 1, // This can be a reactive property based on the current route if needed
      experiment: {
        experimentId: null as string | null,
        metadata: null
      },
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
      elements: [],
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
              "type": "component",
              "inputs": ["pdb_file", "smiles_string"],
              "outputs": ["docked_pdb_file"],
              "description": "Performs docking of protein on ligand."
            },
            {
              "id": "4",
              "name": "GenerateProteins",
              "type": "component",
              "inputs": [],
              "outputs": ["generated_fasta_file"],
              "description": "Generates new protein sequences."
            },
            {
              "id": "5",
              "name": "RunFolding",
              "type": "component",
              "inputs": ["generated_fasta_file"],
              "outputs": ["generated_pdb_file"],
              "description": "Folds generated protein sequences into 3D structures."
            },
            {
              "id": "6",
              "name": "DockingProteinOnProtein",
              "type": "component",
              "inputs": ["pdb_file", "generated_pdb_file"],
              "outputs": ["protein_complex_pdb_file"],
              "description": "Performs docking of generated protein on target protein."
            },
            {
              "id": "7",
              "name": "CalculateBlast",
              "type": "component",
              "inputs": ["protein_complex_pdb_file"],
              "outputs": ["similarity_scores"],
              "description": "Calculates similarity scores between protein complexes."
            }
          ],
          "edges": [
            {"from": "1", "to": "3"},
            {"from": "2", "to": "3"},
            {"from": "4", "to": "5"},
            {"from": "5", "to": "6"},
            {"from": "1", "to": "6"},
            {"from": "6", "to": "7"}
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
      const nodes: Node[] = []; // Define nodes as Node array
      const edges: Edge[] = []; // Define edges as Edge array

      // Helper function to get the depth of a node in the graph
      const getNodeDepth = (nodeId: string): number => {
        const incomingEdges = this.jsonData.workflow.edges.filter(edge => edge.to === nodeId);
        if (incomingEdges.length === 0) return 0; // Node has no incoming edges, depth is 0
        const parentDepths = incomingEdges.map(edge => getNodeDepth(edge.from));
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
      const layerWidth = 200; // Width difference between layers
      const layerHeight = 100; // Height difference between nodes within a layer
      Object.values(nodesByDepth).forEach((nodesAtDepth: Node[], depth: number) => {
        const numNodes = nodesAtDepth.length;
        const startY = 100 + (numNodes * layerHeight) / 2;

        nodesAtDepth.forEach((node, index) => {
          const position = { x: currentX, y: startY - index * layerHeight };
          const nodeData = {
            id: node.id,
            type: node.type,
            label: node.name,
            position,
            data: { description: node.description }
          };

          // Set sourcePosition and targetPosition for special nodes
          if (node.type === 'special') {
            nodeData.sourcePosition = Position.Right;
            nodeData.targetPosition = Position.Left;
            nodeData.specialNodeProps = {
              // Add any special properties needed for the special node component
              // For example:
              // specialProperty: node.specialProperty
            };
          } else {
            // Set default sourcePosition and targetPosition
            nodeData.sourcePosition = Position.Right;
            nodeData.targetPosition = Position.Left;
          }

          nodes.push(nodeData);
        });

        // Update currentX for the next layer
        currentX += layerWidth;
      });

      // Add edges
      this.jsonData.workflow.edges.forEach(edge => {
        edges.push({
          id: `e${edge.from}-${edge.to}`,
          source: edge.from,
          target: edge.to
        });
      });

      this.elements = { nodes, edges };
    },
    addComponent(option: any) {
      // Find the maximum ID among existing nodes
      let maxId = 0;
      this.elements?.nodes.forEach(node => {
        const idNumber = parseInt(node.id);
        if (!isNaN(idNumber) && idNumber > maxId) {
          maxId = idNumber;
        }
      });

      // Increment the maximum ID by one to generate a new unique ID
      const newNodeId = `${maxId + 1}`;

      // Create the new node
      const newNode: Node = {
        id: newNodeId,
        type: option.value,
        label: option.label,
        position: { x: 100, y: 100 },
        data: { description: option.description },
        sourcePosition: Position.Right,
        targetPosition: Position.Left
      };

      // Add the new node to the elements
      this.elements?.nodes.push(newNode);
    },
    onConnect(params) {
      const newEdge = {
        id: `e${params.source}-${params.target}`,
        source: params.source,
        target: params.target,
      };
      this.elements.edges.push(newEdge);
    },
    onDeleteNode(nodeId: String) {
      console.log(nodeId);
      const index = this.elements.nodes.findIndex(node => node.id === nodeId);
      if (index !== -1) {
        this.elements.nodes.splice(index, 1);
      }
    },
    onNodeDragStopHandler(event: any) {
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
    }
  }
}
)
</script>

<style scoped>

body {
  overflow-x: hidden;
}

.map-container {
  width: 800px;
  height: 800px;
}

.workflow {
  height: 100%;
  width: 100%;
}

</style>
