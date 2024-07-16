import { defineStore } from 'pinia';
import {
  InputPropertyErrorResponse,
  JobErrorResponse,
  JobsACommonControllerForJobsManagementService,
  LigandMetadataResponse,
  ProteinMetadataResponse,
} from 'src/refinedApi/client';
import {
  deleteProtein,
  getAllProteinsMetadata,
  uploadProtein,
  updateProteinName,
  uploadLigand,
  deleteLigand,
  getAllLigandsMetadata,
  getComponentState
} from 'src/features/workflow/refinedApi';
import { Edge, Node as FlowNode } from '@vue-flow/core';
import { v4 as uuidv4 } from 'uuid';
import { Notify } from 'quasar';
import {
  getWorkflow,
  sendWorkflowUpdate,
  getExperimentsApi,
  createExperimentApi,
  deleteExperimentApi
} from 'src/features/workflow/refinedApi';
import {
  ComponentModel_Output,
  WorkflowSchemaModel_Input,
  ComponentModel_Input,
  PropertyModel_Output,
  WorkflowComponentModel,
  MappingModel,
  DefaultWorkflowComponentModelValue
} from 'src/refinedApi/client';

// Define custom Node type
export interface Node extends FlowNode {
  id: string;
  name: string;
  type: string;
  inputs: string[];
  outputs: string[];
  description: string;
  jobIds: string[];
  jobs_errors: JobErrorResponse[];
  input_property_errors: InputPropertyErrorResponse[];
  last_exceptions: string[];
  error: string;
  defaults: Array<DefaultWorkflowComponentModelValue>;
}

export const useWorkflowStore = defineStore('workflowStore', {
  state: () => ({
    proteins: [] as ProteinMetadataResponse[],
    ligands: [] as LigandMetadataResponse[],
    elements: {
      nodes: [] as Node[],
      edges: [] as Edge[]
    },
    jobIdsToUpdate: [] as string[],
    workflowId: "" as string,
    allowedTypes: ["Proteins", "Ligands", "DNA"],
    componentOptions: [] as Array<{
      name: string;
      type: string;
      inputs: Record<string, PropertyModel_Output>,
      outputs: Record<string, PropertyModel_Output>,
      description: string
    }>,
    runningComponentIds: [] as string[]
  }),
  getters: {
    workflowIsRunning(state) {
      return state.runningComponentIds.length > 0;
    }
  },
  actions: {
    async createExperiment() {
      return await createExperimentApi();
    },
    async deleteExperiment(experimentId: string) {
      return await deleteExperimentApi(experimentId);
    },
    async getExperiments() {
      return await getExperimentsApi();
    },
    async getAllProteins(experimentId: string) {
      const response = await getAllProteinsMetadata(experimentId);
      this.proteins = response;
    },
    addRunningComponentId(componentId: string) {
      if (!this.runningComponentIds.includes(componentId)) {
        this.runningComponentIds.push(componentId);
      }
    },
    removeRunningComponentId(componentId: string) {
      this.runningComponentIds = this.runningComponentIds.filter(id => id !== componentId);
    },
    clearRunningComponentIds() {
      this.runningComponentIds = [];
    },
    async deleteJob(jobId: string) {
      await JobsACommonControllerForJobsManagementService.deleteJobApiV1JobsJodIdDelete(jobId);
    },
    async uploadProteinToExperiment(
      experimentId: string,
      nodeId: string,
      name?: string,
      fastaFile?: Blob,
      pdbFile?: Blob,
      metaData?: Record<string, string>
    ) {
      try {
        let link: string | undefined;
        if (metaData && 'link' in metaData) {
          link = metaData.link;
        }

        const uploadedProtein = await uploadProtein(
          experimentId,
          name,
          fastaFile,
          pdbFile,
          link // Pass the link if it exists
        );

        // Check if the protein already exists in the proteins array
        const proteinExists = this.proteins.some(protein => protein.id === uploadedProtein.id);
        if (!proteinExists) {
          this.proteins.push(uploadedProtein);
        }

        const existingNode = this.getNodeById(nodeId);

        if (existingNode) {
          if (!existingNode.data.defaults || !Array.isArray(existingNode.data.defaults) || existingNode.data.defaults.length < 1) {
            // Initialize defaults if it doesn't exist or is not an array
            existingNode.data.defaults = [{
              target_path: Object.keys(existingNode.data.inputs || {}), value: [uploadedProtein.id]
            }] as DefaultWorkflowComponentModelValue[];
          } else {
            // Find the default entry, if it exists
            const defaultEntry = existingNode.data.defaults.find(entry =>
              Array.isArray(entry.target_path) &&
              entry.target_path.every((path, index) => path === Object.keys(existingNode.data.inputs || {})[index])
            );

            if (defaultEntry) {
              // Check if the protein id already exists in the defaultEntry value array
              if (!defaultEntry.value.includes(uploadedProtein.id)) {
                // Push the new id to the existing value array
                defaultEntry.value.push(uploadedProtein.id);
              }
            } else {
              // If no matching default entry exists, create a new one
              existingNode.data.defaults.push({
                target_path: Object.keys(existingNode.data.inputs || {}),
                value: [uploadedProtein.id]
              });
            }
          }
          this.sendWorkflowUpdate();
        }
      } catch (error) {
        console.error('Error uploading protein:', error);
      }
    },
    async deleteProteinFromExperiment(nodeId: string, proteinId: string) {
      try {
        await deleteProtein(proteinId);
        this.proteins = this.proteins.filter(protein => protein.id !== proteinId);

        const existingNode = this.getNodeById(nodeId);

        if (existingNode && existingNode.data.defaults && Array.isArray(existingNode.data.defaults)) {
          const defaultEntry = existingNode.data.defaults.find(entry =>
            Array.isArray(entry.target_path) &&
            entry.target_path.every((path, index) => path === Object.keys(existingNode.data.inputs || {})[index])
          );

          if (defaultEntry) {
            defaultEntry.value = defaultEntry.value.filter(id => id !== proteinId);
          }
        }
        this.sendWorkflowUpdate();
      } catch (error) {
        console.error('Error deleting protein:', error);
      }
    },
    async changeProteinName(proteinId: string, newName: string) {
      try {
        await updateProteinName(proteinId, newName);
        const protein = this.proteins.find(protein => protein.id === proteinId);
        if (protein) {
          protein.name = newName;
        }
      } catch (error) {
        console.error('Error deleting protein:', error);
      }
    },
    async getAllLigands(experimentId: string) {
      const response = await getAllLigandsMetadata(experimentId);
      this.ligands = response;
    },
    async uploadLigandToExperiment(
      experimentId: string,
      nodeId: string,
      name?: string,
      smiles?: Blob,
      sdf?: Blob,
      metaData?: Record<string, string>
    ) {
      try {
        const uploadedLigand = await uploadLigand(
          experimentId,
          name,
          smiles,
          sdf
        );

        // Check if the ligand already exists in the ligands array
        const ligandExists = this.ligands.some(ligand => ligand.id === uploadedLigand.id);
        if (!ligandExists) {
          this.ligands.push(uploadedLigand);
        }

        const existingNode = this.getNodeById(nodeId);

        if (existingNode) {
          if (!existingNode.data.defaults) {
            existingNode.data.defaults = [];
          }

          let defaultEntry = existingNode.data.defaults.find(entry =>
            Array.isArray(entry.target_path) &&
            entry.target_path.every((path, index) => path === Object.keys(existingNode.data.inputs || {})[index])
          );

          if (!defaultEntry) {
            defaultEntry = {
              target_path: Object.keys(existingNode.data.inputs || {}),
              value: []
            } as DefaultWorkflowComponentModelValue;
            existingNode.data.defaults.push(defaultEntry);
          }

          // Check if the ligand id already exists in the defaultEntry value array
          if (!defaultEntry.value.includes(uploadedLigand.id)) {
            defaultEntry.value.push(uploadedLigand.id);
          }
        }
        this.sendWorkflowUpdate();
      } catch (error) {
        console.error('Error uploading ligand:', error);
      }
    },
    async deleteLigandFromExperiment(nodeId: string, ligandId: string) {
      try {
        await deleteLigand(ligandId);
        this.ligands = this.ligands.filter(ligand => ligand.id !== ligandId);

        const existingNode = this.getNodeById(nodeId);

        if (existingNode && existingNode.data.defaults && Array.isArray(existingNode.data.defaults)) {
          const defaultEntry = existingNode.data.defaults.find(entry =>
            Array.isArray(entry.target_path) &&
            entry.target_path.every((path, index) => path === Object.keys(existingNode.data.inputs || {})[index])
          );

          if (defaultEntry) {
            defaultEntry.value = defaultEntry.value.filter(id => id !== ligandId);
          }
        }
        this.sendWorkflowUpdate();
      } catch (error) {
        console.error('Error deleting ligand:', error);
      }
    },
    async adjustWorkflow(data: any) {
      // Retain existing nodes
      const existingNodeIds = this.elements.nodes.map(node => node.id);
    
      // Store old to new ID mapping
      const idMapping = new Map<string, string>();
    
      // Process new nodes and generate UUIDs
      const newNodes = data.workflow_components.map(component => {
        if (existingNodeIds.includes(component.id)) {
          return component; // Retain existing component
        } else {
          const newId = uuidv4();
          idMapping.set(component.id, newId); // Map old ID to new ID
          component.id = newId;
    
          // Find component options to get inputs and outputs
          const componentOptions = this.componentOptions.find(opt => opt.name === component.name);
    
          const nodeType = this.allowedTypes.includes(componentOptions?.type || '') ? componentOptions?.type || 'custom' : 'custom';
    
          // Create the new node
          const newNode = {
            id: newId,
            name: component.name,
            description: component.description,
            type: nodeType,
            data: {
              description: componentOptions ? componentOptions.description : '',
              inputs: componentOptions ? componentOptions.inputs : {},
              outputs: componentOptions ? componentOptions.outputs : {},
              jobIds: [],
              jobs_errors: [],
              input_property_errors: [],
              last_exceptions: [],
              draggable: false,
              defaults: [],
              error: null
            },
            position: {
              x: component.x !== undefined ? component.x : 100,
              y: component.y !== undefined ? component.y : 100
            }
          };
    
          this.elements.nodes.push(newNode);
    
          return component;
        }
      });
    
      // Adjust connections (edges)
      const newEdges = data.workflow_components.flatMap(component =>
        component.connections.map(conn => {
          // Adjust source_component_id if it matches an old ID
          const adjustedSourceId = idMapping.get(conn.source_component_id) || conn.source_component_id;
    
          return {
            id: `e${adjustedSourceId}-to-${component.id}`,
            source: adjustedSourceId,
            target: component.id,
            sourceHandle: `${adjustedSourceId}-output-${conn.source_output}`,
            targetHandle: `${component.id}-input-${conn.target_input}`,
            type: 'custom'
          };
        })
      );
    
      // Remove duplicate edges by using a Set
      const uniqueEdges = Array.from(new Set(newEdges.map(edge => edge.id)))
        .map(id => newEdges.find(edge => edge.id === id));
    
      // Remove edges that are no longer present
      this.elements.edges = this.elements.edges.filter(edge => uniqueEdges.some(newEdge => newEdge.id === edge.id));
      
      // Add new edges
      uniqueEdges.forEach(edge => {
        if (!this.elements.edges.some(existingEdge => existingEdge.id === edge.id)) {
          this.elements.edges.push(edge);
        }
      });
    
      // Remove nodes that are no longer present
      const newNodeIds = newNodes.map(node => node.id);
      this.elements.nodes = this.elements.nodes.filter(node => newNodeIds.includes(node.id));
    
      // Update the state
      await this.sendWorkflowUpdate();

      await this.fetchWorkflow(this.workflowId);
    },
    // Helper function to fetch component options
    getComponentOptionsByName(name: string) {
      return this.componentOptions.find(opt => opt.name === name) || {};
    },
    async fetchWorkflow(workflowId: string) {
      this.workflowId = workflowId;
      try {
        const workflow = await getWorkflow(workflowId);
        if (!workflow) {
          console.error("Failed to load workflow");
          return;
        }

        const nodes: Node[] = [];
        const edges: Edge[] = [];

        const componentIOMap: {
          [key: string]: {
            inputs: Record<string, PropertyModel_Output>,
            outputs: Record<string, PropertyModel_Output>,
            type: string,
            description: string
          }
        } = {};
        workflow.components.forEach(component => {
          const inputs = component.input;
          const outputs = component.output;
          const description = this.getDescriptionString(component);
          componentIOMap[component.name] = { inputs, outputs, type: component.name, description };
        });

        for (const component of workflow.workflow_components) {
          const { inputs, outputs, type, description } = componentIOMap[component.name] || {
            inputs: [],
            outputs: [],
            type: '',
            description: ''
          };
          const nodeType = this.allowedTypes.includes(type) ? type : "custom";

          const componentState = await getComponentState(component.component_id);

          const nodeData: Node = {
            id: component.component_id,
            name: component.name,
            description: '',
            type: nodeType,
            jobIds: componentState.job_ids,
            jobs_errors: componentState.jobs_errors,
            input_property_errors: componentState.input_property_errors,
            last_exceptions: componentState.last_exceptions,
            data: {
              description: description,
              inputs,
              outputs,
              draggable: false,
              defaults: component.defaults,
              error: component.error
            },
            position: {
              x: component.x !== undefined ? component.x : 100,
              y: component.y !== undefined ? component.y : 100
            }
          };
          nodes.push(nodeData);
        }

        workflow.workflow_components.forEach(component => {
          component.mappings?.forEach(mapping => {
            edges.push({
              id: `e${mapping.source_component_id}-to-${component.component_id}`,
              source: mapping.source_component_id,
              target: component.component_id,
              sourceHandle: `${mapping.source_component_id}-output-${mapping.source_path[0]}`,
              targetHandle: `${component.component_id}-input-${mapping.target_path[0]}`,
              type: "custom",
              data: { text: mapping.error }
            });
          });
        });

        this.elements.nodes = nodes;
        this.elements.edges = edges;

        this.componentOptions = workflow.components.map(component => ({
          name: component.name,
          type: component.name,
          inputs: component.input,
          outputs: component.output,
          description: component.description
        }));
      } catch (error) {
        console.error("Error fetching workflow:", error);
      }
    },
    getDescriptionString(component: ComponentModel_Output): string {
      const inputKeys = Object.keys(component.input || {});
      const description = inputKeys.length > 0 ? component.input[inputKeys[0]].description : null;
      return description || "No description available";
    },
    async sendWorkflowUpdate() {
      const workflowUpdate: WorkflowSchemaModel_Input = {
        workflow_id: this.workflowId,
        components: this.componentOptions.map(option => ({
          name: option.name,
          input: option.inputs,
          output: option.outputs
        }) as ComponentModel_Input),
        workflow_components: this.elements.nodes.map(node => ({
          name: node.name,
          component_id: node.id,
          error: node.data.error,
          x: node.position.x,
          y: node.position.y,
          mappings: this.elements.edges
            .filter(edge => edge.target === node.id)
            .map(edge => ({
              source_path: [edge.sourceHandle?.split('-output-')[1]],
              target_path: [edge.targetHandle?.split('-input-')[1]],
              source_component_id: edge.sourceHandle?.split('-output-')[0],
              error: edge.data?.text
            }) as MappingModel),
          defaults: node.data.defaults,
        }) as WorkflowComponentModel),
        error: null,
        valid: true
      };
      try {
        await sendWorkflowUpdate(workflowUpdate);
        console.log('Workflow updated successfully');
      } catch (error) {
        console.error('Error updating workflow:', error);
      }
    },
    addComponent(option: any) {
      // Generate a new UUID for the node ID
      const newNodeId = uuidv4();
      const nodeType = this.allowedTypes.includes(option.type) ? option.type : "custom";

      // Create the new node
      const newNode: Node = {
        id: newNodeId,
        name: option.name,
        type: nodeType,
        data: {
          description: option.description,
          inputs: option.inputs,
          outputs: option.outputs,
          jobIds: [],
          jobs_errors: [],
          input_property_errors: [],
          last_exceptions: [],
          draggable: false,
          defaults: [],
          error: null
        },
        position: { x: 100, y: 100 }, // Default position
      };

      // Add the new node to the elements
      this.elements.nodes.push(newNode);

      // Send the workflow update
      this.sendWorkflowUpdate();
    },
    onConnect(params: { source: string, target: string, sourceHandle: string, targetHandle: string }) {
      // Extract node IDs and handle types
      const sourceNodeId = params.sourceHandle.split('-output-')[0];
      const targetNodeId = params.targetHandle.split('-input-')[0];
      const sourceHandleType = params.sourceHandle.includes('-output-') ? 'output' : 'input';
      const targetHandleType = params.targetHandle.includes('-input-') ? 'input' : 'output';

      // Validation 1: Only outputs and inputs can be connected (and not vice versa)
      if (sourceHandleType !== 'output' || targetHandleType !== 'input') {
        Notify.create({
          type: 'negative',
          message: 'Only outputs and inputs can be connected.',
        });
        return;
      }

      // Validation 2: Outputs and inputs of the same node can't be connected
      if (sourceNodeId === targetNodeId) {
        Notify.create({
          type: 'negative',
          message: 'Cannot connect outputs and inputs of the same node.',
        });
        return;
      }

      // Validation 3: Check if the same edge already exists
      const edgeExists = this.elements.edges.some(edge =>
        edge.source === sourceNodeId &&
        edge.target === targetNodeId &&
        edge.sourceHandle === params.sourceHandle &&
        edge.targetHandle === params.targetHandle
      );

      if (edgeExists) {
        Notify.create({
          type: 'negative',
          message: 'This connection already exists.',
        });
        return;
      }

      // If all validations pass, create and add the new edge
      const newEdge: Edge = {
        id: `e${params.sourceHandle}-to-${params.targetHandle}`,
        source: sourceNodeId,
        target: targetNodeId,
        sourceHandle: params.sourceHandle,
        targetHandle: params.targetHandle,
        type: "custom"
      };

      this.elements.edges.push(newEdge);

      // Send the workflow update
      this.sendWorkflowUpdate();
    },
    onEdgeRemove(edgeId: string) {
      const index = this.elements.edges.findIndex(edge => edge.id === edgeId);
      if (index !== -1) {
        this.elements.edges.splice(index, 1);

        this.sendWorkflowUpdate();
      }
    },
    onDeleteNode(nodeId: string) {
      const index = this.elements.nodes.findIndex(node => node.id === nodeId);
      if (index !== -1) {
        this.elements.nodes.splice(index, 1);
        this.elements.edges = this.elements.edges.filter(edge => edge.source !== nodeId && edge.target !== nodeId);

        // Send the workflow update
        this.sendWorkflowUpdate();
      }
    },
    onNodeDragStopHandler(event: { node: Node }) {
      const id = event.node.id;
      const newPosition = event.node.position;

      // Find the node with the given id and update its position
      const nodeToUpdate = this.elements.nodes.find((node: Node) => node.id === id);
      if (nodeToUpdate) {
        nodeToUpdate.position = newPosition;

        this.sendWorkflowUpdate();
      }
    },
    getNodeById(nodeId: string) {
      const node = this.elements.nodes.find(node => node.id === nodeId) || null;
      return node;
    },
    addJobIdToUpdate(jobId: string) {
      if (!this.jobIdsToUpdate.includes(jobId)) {
        this.jobIdsToUpdate.push(jobId);
      }
    },
    removeJobIdToUpdate(jobId: string) {
      this.jobIdsToUpdate = this.jobIdsToUpdate.filter(id => id !== jobId);
    },
    async pollWorkflow() {
      try {
        const workflow = await getWorkflow(this.workflowId);
        if (!workflow) {
          console.error("Failed to load workflow");
          return;
        }

        for (const workflow_component of workflow.workflow_components) {
          const existingNode = this.getNodeById(workflow_component.component_id);
          if (existingNode) {
            const componentState = await getComponentState(workflow_component.component_id);
            existingNode.data.jobIds = componentState.job_ids;
            existingNode.data.jobs_errors = componentState.jobs_errors;
            existingNode.data.input_property_errors = componentState.input_property_errors;
            existingNode.data.last_exceptions = componentState.last_exceptions;
            existingNode.data.error = workflow_component.error;
          }
        }

        const updatedEdges: Edge[] = [];
        workflow.workflow_components.forEach(component => {
          component.mappings?.forEach(mapping => {
            const existingEdge = this.elements.edges.find(edge => edge.source === mapping.source_component_id && edge.target === component.component_id);
            if (existingEdge && existingEdge.data) {
              existingEdge.data.text = mapping.error;
            }
          });
        });

        // Re-assign edges to ensure the VueFlow component updates
        this.elements.edges = [...this.elements.edges.filter(edge => !updatedEdges.some(updatedEdge => updatedEdge.id === edge.id)), ...updatedEdges];
      } catch (error) {
        console.error('Error polling workflow:', error);
      }
    }
  },
});
