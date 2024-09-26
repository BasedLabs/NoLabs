import {defineStore} from 'pinia';
import {
  ComponentSchema,
  ComponentSchemaTemplate_Input,
  ComponentSchemaTemplate_Output, ComponentStateEnum,
  DefaultSchema,
  JobsACommonControllerForJobsManagementService,
  LigandMetadataResponse,
  MappingSchema,
  PropertyErrorResponse,
  PropertySchema_Input,
  PropertySchema_Output,
  ProteinMetadataResponse,
  WorkflowSchema_Input,
} from 'src/refinedApi/client';
import {
  createExperimentApi,
  deleteExperimentApi,
  deleteLigand,
  deleteProtein,
  getAllLigandsMetadata,
  getAllProteinsMetadata,
  getComponentState,
  getExperimentsApi,
  getWorkflow,
  sendWorkflowUpdate,
  updateProteinName,
  uploadLigand,
  uploadProtein
} from 'src/features/workflow/refinedApi';
import {Node as FlowNode} from '@vue-flow/core';
import {v4 as uuidv4} from 'uuid';
import {Notify} from 'quasar';
import {BiobuddyWorkflowAdjustmentData} from "../../biobuddy/biobuddyTypes";

// Define a new NodeData type
interface NodeData {
  description: string;
  inputs: Record<string, PropertySchema_Input>;
  outputs: Record<string, PropertySchema_Output>;
  jobIds: string[];
  jobIdsToUpdate?: string[];
  input_property_errors: PropertyErrorResponse[];
  state: ComponentStateEnum | null,
  stateMessage: string;
  draggable: boolean;
  defaults?: DefaultSchema[];
  error?: string | null;
}

// Update the Node interface to use the new NodeData type
export interface Node extends FlowNode {
  id: string;
  name: string;
  type: string;
  data: NodeData; // Use the new NodeData type here
}

export interface Edge {
  id: string;
  source: string;
  target: string;
  sourceHandle: string;
  targetHandle: string;
  type: string;
  data?: { text: string | null | undefined; };
}

export const useWorkflowStore = defineStore('workflowStore', {
  state: () => ({
    proteins: [] as ProteinMetadataResponse[],
    ligands: [] as LigandMetadataResponse[],
    elements: {
      nodes: [] as Node[],
      edges: [] as Edge[]
    },
    experimentId: "" as string,
    allowedTypes: ["Proteins", "Ligands", "DNA"],
    componentOptions: [] as Array<{
      name: string;
      type: string;
      inputs: Record<string, PropertySchema_Input>,
      outputs: Record<string, PropertySchema_Output>,
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
      this.proteins = await getAllProteinsMetadata(experimentId);
    },
    async deleteJob(jobId: string, componentId: string) {
      try {
        // Delete the job via the API
        await JobsACommonControllerForJobsManagementService.deleteJobApiV1JobsJodIdDelete(jobId);
        console.log(`Job ${jobId} deleted successfully`);

        // Find the component by componentId
        const component = this.getNodeById(componentId);

        if (component && component.data && component.data.jobIds) {
          // Remove the jobId from component.data.jobIds
          component.data.jobIds = component.data.jobIds.filter((id: string) => id !== jobId);
          console.log(`Job ${jobId} removed from component ${componentId}`);
        } else {
          console.error(`Component with ID ${componentId} not found or has no jobIds`);
        }
      } catch (error) {
        console.error(`Error deleting job ${jobId}:`, error);
      }
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
            }] as Array<DefaultSchema>;
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
    setInputValue(nodeId: string, inputName: string, inputValue: any) {
      const existingNode = this.getNodeById(nodeId);
      if (existingNode) {
        if (!existingNode.data.defaults || !Array.isArray(existingNode.data.defaults) || existingNode.data.defaults.length < 1) {
          // Initialize defaults if it doesn't exist or is not an array
          existingNode.data.defaults = [{
            target_path: [inputName],
            value: inputValue
          }];
        } else {
          // Check if the inputName already exists in the defaults array
          const existingDefault = existingNode.data.defaults.find(def => def.target_path.includes(inputName));
          if (existingDefault) {
            existingDefault.value = inputValue;
          } else {
            // Add a new default entry if it doesn't exist
            existingNode.data.defaults.push({
              target_path: [inputName],
              value: inputValue
            });
          }
        }
        this.sendWorkflowUpdate();
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
            defaultEntry.value = defaultEntry.value.filter((id: string) => id !== proteinId);
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
      this.ligands = await getAllLigandsMetadata(experimentId);
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

          // Find or create the default entry
          let defaultEntry = existingNode.data.defaults.find(entry =>
            Array.isArray(entry.target_path) &&
            entry.target_path.every((path, index) => path === Object.keys(existingNode.data.inputs || {})[index])
          );

          if (!defaultEntry) {
            defaultEntry = {
              target_path: Object.keys(existingNode.data.inputs || {}),
              value: []
            } as DefaultSchema;
            existingNode.data.defaults.push(defaultEntry);
          }

          // Initialize the value array if it's not an array
          if (!Array.isArray(defaultEntry.value)) {
            defaultEntry.value = [];
          }

          // Check if the ligand id already exists in the defaultEntry value array
          if (!defaultEntry.value.includes(uploadedLigand.id)) {
            // Push the new id to the existing value array
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
            defaultEntry.value = defaultEntry?.value.filter((id: string) => id !== ligandId);
          }
        }
        this.sendWorkflowUpdate();
      } catch (error) {
        console.error('Error deleting ligand:', error);
      }
    },
    updateJobStatus(componentId: string, jobId: string) {
      const component = this.getNodeById(componentId);

      if (component && component.data && component.data.jobIdsToUpdate) {
        // Find the job in the component's jobIds and update its status
        if (!component.data.jobIdsToUpdate.includes(jobId)) {
          component.data.jobIdsToUpdate.push(jobId);
        }
      } else {
        console.error(`Component ${componentId} not found or has no jobIds`);
      }
    },
    async adjustWorkflow(data: BiobuddyWorkflowAdjustmentData) {
      // Retain existing nodes
      const existingNodeIds = this.elements.nodes.map(node => node.id);

      // Store old to new ID mapping
      const idMapping = new Map<string, string>();

      // Process new nodes and generate UUIDs
      const newNodes = data.components.map(component => {
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
              jobIdsToUpdate: [],
              jobs_errors: [],
              input_property_errors: [],
              state: null,
              stateMessage: '',
              draggable: false,
              defaults: [],
              error: null
            } as NodeData,
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
      const newEdges = data.components.flatMap(component =>
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
      this.elements.edges = this.elements.edges.filter(edge => uniqueEdges.some(newEdge => newEdge?.id === edge.id));

      // Add new edges
      uniqueEdges.forEach(edge => {
        if (edge && !this.elements.edges.some(existingEdge => existingEdge.id === edge?.id)) {
          this.elements.edges.push(edge);
        }
      });

      // Remove nodes that are no longer present
      const newNodeIds = newNodes.map(node => node.id);
      this.elements.nodes = this.elements.nodes.filter(node => newNodeIds.includes(node.id));

      // Update the state
      await this.sendWorkflowUpdate();

      await this.fetchWorkflow(this.experimentId);
    },
    async fetchWorkflow(experimentId: string) {
      this.experimentId = experimentId;
      try {
        const workflow = await getWorkflow(experimentId);
        if (!workflow) {
          console.error("Failed to load workflow");
          return;
        }

        const nodes: Node[] = [];
        const edges: Edge[] = [];

        const componentIOMap: {
          [key: string]: {
            inputs: Record<string, PropertySchema_Output>,
            outputs: Record<string, PropertySchema_Output>,
            type: string,
            description: string
          }
        } = {};
        workflow.component_templates.forEach(component => {
          const inputs = component.input;
          const outputs = component.output;
          const description = this.getDescriptionString(component);
          componentIOMap[component.name] = { inputs, outputs, type: component.name, description };
        });

        for (const component of workflow.components) {
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
            type: nodeType,
            data: {
              description: description,
              inputs,
              outputs,
              draggable: false,
              defaults: component.defaults,
              error: component.error,
              jobIds: componentState.job_ids,
              jobIdsToUpdate: [],
              input_property_errors: componentState.input_errors,
              state: componentState.state,
              stateMessage: componentState.state_message,
            } as NodeData,
            position: {
              x: component.x !== undefined ? component.x : 100,
              y: component.y !== undefined ? component.y : 100
            }
          };
          nodes.push(nodeData);
        }

        workflow.components.forEach(component => {
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

        this.componentOptions = workflow.component_templates.map(component => ({
          name: component.name,
          type: component.name,
          inputs: component.input,
          outputs: component.output,
          description: component.description as string
        }));
      } catch (error) {
        console.error("Error fetching workflow:", error);
      }
    },
    getDescriptionString(component: ComponentSchemaTemplate_Output): string {
      const inputKeys = Object.keys(component.input || {});
      const description = inputKeys.length > 0 ? component.input[inputKeys[0]].description : null;
      return description || "No description available";
    },
    async sendWorkflowUpdate() {
      const workflowUpdate: WorkflowSchema_Input = {
        experiment_id: this.experimentId,
        component_templates: this.componentOptions.map(option => ({
          name: option.name,
          input: option.inputs,
          output: option.outputs
        }) as ComponentSchemaTemplate_Input),
        components: this.elements.nodes.map(node => ({
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
            }) as MappingSchema),
          defaults: node.data.defaults,
        }) as ComponentSchema),
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
          jobIdsToUpdate: [],
          input_property_errors: [],
          state: null,
          stateMessage: '',
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
    updateComponentState(componentId: string, newState: ComponentStateEnum) {
      const component = this.getNodeById(componentId);  // Assuming you have getNodeById method to find the component

      if (component) {
        // Update the state of the component
        component.data.state = newState;

        // Manage the runningComponentIds list based on the new state
        if (newState === ComponentStateEnum.RUNNING) {
          // Add the component to the list of running components if it's not already there
          if (!this.runningComponentIds.includes(componentId)) {
            this.runningComponentIds.push(componentId);
            console.log(`Component ${componentId} added to running components.`);
          }
        } else if (newState === ComponentStateEnum.COMPLETED) {
          // Remove the component from the running list once it has completed
          this.runningComponentIds = this.runningComponentIds.filter(id => id !== componentId);
          console.log(`Component ${componentId} removed from running components.`);
        }

        console.log(`Component ${componentId} updated to state: ${newState}`);
      } else {
        console.error(`Component with ID ${componentId} not found`);
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
    getNodeById(nodeId: string): Node | null {
      return this.elements.nodes.find(node => node.id === nodeId) || null;
    },
    removeJobIdToUpdate(componentId: string, jobId: string) {
      const component = this.getNodeById(componentId);
      if (component && component.data && component.data.jobIdsToUpdate){
        component.data.jobIdsToUpdate = component.data.jobIdsToUpdate.filter(id => id !== jobId);
      }
    },
    async adjustComponentJobsList(componentId: string, jobIds: string[]) {
      try {
        // Find the node by componentId
        const existingNode = this.getNodeById(componentId);

        if (!existingNode) {
          console.error(`Component with ID ${componentId} not found`);
          return;
        }

        // Update the jobIds of the existing node's data
        existingNode.data.jobIds = jobIds;

        // Optionally, you can trigger any reactivity updates or other logic here
        console.log(`Updated job list for component ${componentId}:`, jobIds);

      } catch (error) {
        console.error(`Error adjusting component jobs list for component ${componentId}:`, error);
      }
    }

  },
});
