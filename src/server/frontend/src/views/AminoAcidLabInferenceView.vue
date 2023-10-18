<script>
export default {
    data() {
        return {
            tabId: 'localisation'
        }
    },
    methods: {
        onTabClick (tabId) {
            this.tabId = tabId;

            if (this.tabId === 'protein3dViewer') {
                this.$refs.proteinViewerInstance.value.render();
            }
        }
    }
}
</script>

<template>
    <nav>
        <div class="nav nav-tabs justify-content-center" id="nav-tab" role="tablist">
            <button class="nav-link" :class="tabId === 'localisation' ? 'active' : ''" type="button"
                @click="onTabClick('localisation')" role="tab" :aria-selected="tabId === 'localisation'">
                Localisation
            </button>
            <button class="nav-link" :class="tabId === 'protein3dViewer' ? 'active' : ''"
                @click="onTabClick('protein3dViewer')" type="button" role="tab"
                :aria-selected="tabId === 'protein3dViewer'">
                Protein 3D viewer
            </button>
            <button class="nav-link" :class="tabId === 'geneOntology' ? 'active' : ''"
                @click="onTabClick('geneOntology')" type="button" role="tab"
                :aria-selected="tabId === 'geneOntology'">
                Gene ontology
            </button>
            <button class="nav-link" :class="tabId === 'solubility' ? 'active' : ''"
                @click="onTabClick('solubility')" type="button" role="tab"
                :aria-selected="tabId === 'solubility'">
                Solubility
            </button>
        </div>
    </nav>
    <div class="tab-content" id="nav-tabContent">
        <div class="tab-pane fade mt-1" :class="tabId === 'localisation' ? 'show active' : ''" role="tabpanel">
            <ProteinLocalisation />
        </div>
        <div class="tab-pane fade mt-1" :class="tabId === 'protein3dViewer' ? 'show active' : ''" role="tabpanel">
            <ProteinViewer ref="proteinViewerInstance" />
        </div>
        <div class="tab-pane fade mt-1" :class="tabId === 'geneOntology' ? 'show active' : ''" role="tabpanel">
            <GeneOntology />
        </div>
        <div class="tab-pane fade mt-1" :class="tabId === 'solubility' ? 'show active' : ''" role="tabpanel">
            <ProteinSolubility />
        </div>
    </div>
</template>