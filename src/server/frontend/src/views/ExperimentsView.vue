<script>
import { io } from "socket.io-client";
import { h } from "vue";
const socket = io("http://127.0.0.1:5000");

export default {
    props: ['state', 'api'],
    data() {
        return {
            serverLogs: 'Loading . . .',
            selectedExperiment: null,
            showSideMenu: false,
            renderSideMenu: false
        }
    },
    methods: {
        startLoader() {
            return this.$loading.show(null, {default: () => h('div', {}, this.serverLogs )});
        },
        async loadExperiments() {
            const loader = this.startLoader();
            await this.api.getAllExperiments();
            loader.hide();
        },
        async selectExperiment(experiment) {
            const loader = this.startLoader();
            await this.api.loadExperiment(experiment);
            loader.hide();
            this.selectedExperiment = experiment;
        },
        async addExperiment() {
            const loader = this.startLoader();
            await this.api.addExperiment();
            loader.hide();
        },
        async deleteExperiment(experiment) {
            const loader = this.startLoader();
            await this.api.deleteExperiment(experiment);
            loader.hide();
        },
        async onFormSubmit(data) {
            const loader = this.startLoader();
            await this.api.inference({ form: data.target, experiment: this.state.experiment });
            loader.hide();
        },
        isCurrentExperiment(experiment) {
            return this.state.experiment && (this.state.experiment.id == experiment.id)
        },
        async changeExperimentName(evt, experiment) {
            experiment.name = evt.target.innerText;
            await this.api.changeExperimentName(experiment);
        },
        toggleSideMenu() {
            if (this.showSideMenu) {
                // Start the hide transition
                this.showSideMenu = false;
                setTimeout(() => {
                this.renderSideMenu = false;
                }, 300); // Delay should match the CSS transition time
            } else {
                // Start the show transition
                this.renderSideMenu = true;
                this.$nextTick(() => {
                this.showSideMenu = true;
                });
            }
        },
        finishEditing(event) {
            // Remove focus from the element when Enter is pressed
            event.target.blur();
        },
    },
    async mounted() {
        await this.loadExperiments();
        socket.on("get-logs", (...args) => {
            this.serverLogs = args[0].response;
        });
    },
    computed: {
        experimentEmpty() {
            if (!this.state.experiment || !this.state.experiment.data)
                return true;

            if (typeof this.state.experiment.data === 'Array')
                return this.state.experiment.data.length === 0;

            return Object.keys(this.state.experiment.data).length === 0;
        },
        experimentSelected() {
            return !!this.state.experiment;
        }
    }
}
</script>

<template>
    <!-- Center Container -->
    <div class="experiments-menu">
       <div class="side-menu-container">
           <div v-if="selectedExperiment" class="burger-icon" @click="toggleSideMenu">
                   &#9776; <!-- Representing the burger icon -->
           </div>
           <div v-if="renderSideMenu" class="side-menu">
               <div v-for="experiment in state.experiments" :key="experiment.id" class="experiment" @click="selectExperiment(experiment)">
                <h3 
                    contenteditable="true" 
                    @blur="changeExperimentName($event, experiment)"
                    @keyup.enter="finishEditing($event)"
                >{{ experiment.name }}</h3>
               <p>Last Modified: {{ experiment.date }}</p>
               <button class="delete-button" @click.stop="deleteExperiment(experiment)">Delete</button>
               </div>
           </div>
       </div>
       <div v-if="!selectedExperiment" class="main-content">
            <slot name="labTitle"></slot>
            <div v-if="!selectedExperiment" class="add-experiments">
               <button type="button" class="btn btn-primary btn-md" @click.stop="addExperiment()">Add Experiment</button>
            </div>
            <div class="container-fluid row experiments-container">
               <div class="col-md-4"></div>
               <div v-for="experiment in state.experiments" :key="experiment.id" class="col-md-4 experiment" @click="selectExperiment(experiment)">
                    <h3>{{ experiment.name }}</h3>
                    <p>Last Modified: {{ experiment.date }}</p>
                    <div class="tags">
                            <span v-for="type in ['gene ontology', 'folding', 'solubility']" :key="type" class="tag"> {{ type }} </span>
                    </div>
                </div>
                <div class="col-md-4"></div>
            </div>
        </div>
        <div class="col-md-8" v-if="selectedExperiment">
            <div class="text-center m-5">
                <slot name="labTitle"></slot>
                <div class="row" v-if="selectedExperiment">
                    <div class="col-md-12">
                        <slot name="labForm" :onFormSubmit="onFormSubmit"></slot>
                    </div>
                </div>
                <div class="row" v-else>
                    <div class="col-md-12">
                        <hr />
                        <h4>Add or select experiment</h4>
                    </div>
                </div>
                <div id="resultContainer" v-if="!experimentEmpty">
                    <div class="tab-content" id="nav-tabContent">
                        <div class="tab-pane fade mt-1 show active" role="tabpanel">
                            <slot name="lab" :key="state.experiment.id" :experiment="state.experiment"></slot>
                        </div>
                    </div>
                </div>
            </div>
        </div>
   </div>
</template>

<style>
.experiments-menu {
  display: flex;
  justify-content: left;
  align-items: left;
  align-content: center;
  height: 100vh;
  width: 100vw;
}

.add-experiments {
  text-align: center; /* Center the button */
  margin-top: 20px; /* Add some space above the button */
}

.experiments-container {
  margin-left: 250px;
  min-width: 40vw;
  max-width: 80vw;
  max-height: 500px; /* Set a maximum height */
  margin: 20px;
  overflow-y: auto; /* Enables vertical scrolling */
  border: 1px solid #ccc; /* Optional: adds a border around the container */
  border-radius: 5px; /* Optional: rounds the corners */
  padding: 10px; /* Optional: adds some padding inside the container */
}

.experiment {
  background-color: #f0f0f0;
  margin: 10px 0;
  padding: 10px;
  border: 1px solid #ddd;
  border-radius: 5px;
  cursor: pointer;
}

.delete-button {
  cursor: pointer;
  background-color: rgb(232, 59, 59);
  color: #f0f0f0;
  border: none;
  padding: 5px 10px;
  border-radius: 10px;
  /* Additional styling as needed */
}


.tag {
  display: inline-block;
  background-color: #007bff; /* Example background color */
  color: white;
  padding: 5px 10px;
  margin-right: 5px;
  border-radius: 15px; /* Creates the pill shape */
  font-size: 0.8em;
}

.side-menu-container {
  display: flex;
  justify-content: center;
  align-content: center;
}

.main-content {
  margin-top: 40px;
  display: flex;
  align-items: center;
  flex-direction: column; /* Stacks children vertically */
  width: 100vw; /* Takes the full width available */
}

.burger-icon {
  cursor: pointer;
  position: relative;
  margin-left: 20px;
  right: 10px; /* Adjust as needed */
  top: 10px; /* Adjust as needed */
  z-index: 100; /* Ensures it's above other elements in the side menu */
  /* Additional styling for the burger icon */
}

.side-menu {
  min-width: 40vw;
  max-width: 80vw;
  max-height: 500px; /* Set a maximum height */
  overflow-y: auto; /* Enables vertical scrolling */
  border: 1px solid #ccc; /* Optional: adds a border around the container */
  border-radius: 5px; /* Optional: rounds the corners */
  padding: 10px;
  left: 0;
  top: 0;
  bottom: 0;
  background-color: #322d2d;
  transition: transform 0.3s ease-in-out; /* Animation */
  transform: translateX(-100%); /* Initial state */
}

.experiments-menu:not(.show-menu) .side-menu {
  transform: translateX(0); /* Slide out */
}

.menu-item {
  padding: 10px;
  border-bottom: 1px solid #ccc;
  cursor: pointer;
}

@media (prefers-color-scheme: dark) {
    .experiment {
        background-color: rgb(41, 33, 33);
        margin: 10px 0;
        padding: 10px;
        border: 1px solid #ddd;
        border-radius: 5px;
        cursor: pointer;
    }
}
</style>