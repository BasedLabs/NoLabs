<script setup>
import ProteinLocalisation from '../components/AminoAcidLab/ProteinLocalisation.vue';
import ProteinViewer from '../components/AminoAcidLab/ProteinViewer.vue';
import GeneOntology from '../components/AminoAcidLab/GeneOntology.vue';
import ProteinSolubility from '../components/AminoAcidLab/ProteinSolubility.vue';
import { onMounted, reactive, defineEmits, ref } from 'vue';
import store, {api} from '../storage';

const proteinViewerInstance = ref(null);
let pulledInference = ref(false);

const data = reactive({
    tabId: 'localisation'
});

const onTabClick = (tabId) => {
    data.tabId = tabId;

    if (tabId === 'protein3dViewer') {
        proteinViewerInstance.value.render();
    }
}

const onFormSubmit = async (data) => {
    console.log('pulled')
    await api.aminoAcidLab.inference(data);
    pulledInference.value = true;
}

</script>

<template>
    <div class="text-center container">
        <div class="row" :class="pulledInference ? 'invisible' : ''">
            <div class="col-md-12">
                <form enctype="multipart/form-data" id="inferenceInputForm" v-on:submit.prevent="onFormSubmit">
                    <div class="row justify-content-center">
                        <div class="col-md-6">
                            <label for="inputSequence" class="col-form-label fs-4">Paste amino-acid
                                sequence</label>
                            <input type="text" pattern="[SMIREATZLHJXGWKBDPCFNQYV]*" class="form-control"
                                value="AAACGAGGCAA" id="inputSequence" name="inputSequence"
                                title="Sequence is only allowed to have {'S', 'M', 'I', 'R', 'E', 'A', 'T', 'Z', 'L', 'H', 'J', 'X', 'G', 'W', 'K', 'B', 'D', 'P', 'C', 'F', 'N', 'Q', 'Y', 'V'} tokens.">
                        </div>
                    </div>
                    <div class="row justify-content-center">
                        <div class="col-md-6">
                            <label for="inputSequenceFile" class="col-form-label fs-4">Or paste a .FASTA file
                                with
                                amino-acid
                                sequence</label>
                            <input type="file" multiple="multiple" accept="text/x-fasta" class="form-control"
                                id="inputSequenceFile" name="inputSequenceFile">
                        </div>
                    </div>
                    <div class="row justify-content-center m-4">
                        <div class="col-md-4">
                            <button type="submit" id="submitInference" class="btn btn-lg btn-primary">
                                Submit
                            </button>
                        </div>
                    </div>
                </form>
            </div>
        </div>

        <div class="text-center invisible" id="spinner">
            <p class="inference-components-title">Processing. It can take more than 10 mins depending on your
                PC</p>
            <div class="spinner-border" role="status">
                <span class="visually-hidden">Processing</span>
            </div>
        </div>
        <div class="row" id="resultContainer" v-if="pulledInference">
            <nav>
                <div class="nav nav-tabs justify-content-center" id="nav-tab" role="tablist">
                    <button class="nav-link" :class="data.tabId === 'localisation' ? 'active' : ''" type="button"
                        @click="onTabClick('localisation')" role="tab" :aria-selected="data.tabId === 'localisation'">
                        Localisation
                    </button>
                    <button class="nav-link" :class="data.tabId === 'protein3dViewer' ? 'active' : ''"
                        @click="onTabClick('protein3dViewer')" type="button" role="tab"
                        :aria-selected="data.tabId === 'protein3dViewer'">
                        Protein 3D viewer
                    </button>
                    <button class="nav-link" :class="data.tabId === 'geneOntology' ? 'active' : ''"
                        @click="onTabClick('geneOntology')" type="button" role="tab"
                        :aria-selected="data.tabId === 'geneOntology'">
                        Gene ontology
                    </button>
                    <button class="nav-link" :class="data.tabId === 'solubility' ? 'active' : ''"
                        @click="onTabClick('solubility')" type="button" role="tab"
                        :aria-selected="data.tabId === 'solubility'">
                        Solubility
                    </button>
                </div>
            </nav>
            <div class="tab-content" id="nav-tabContent">
                <div class="tab-pane fade mt-1" :class="data.tabId === 'localisation' ? 'show active' : ''" role="tabpanel">
                    <ProteinLocalisation />
                </div>
                <div class="tab-pane fade mt-1" :class="data.tabId === 'protein3dViewer' ? 'show active' : ''" role="tabpanel">
                    <ProteinViewer ref="proteinViewerInstance" />
                </div>
                <div class="tab-pane fade mt-1" :class="data.tabId === 'geneOntology' ? 'show active' : ''" role="tabpanel">
                    <GeneOntology />
                </div>
                <div class="tab-pane fade mt-1" :class="data.tabId === 'solubility' ? 'show active' : ''" role="tabpanel">
                    <ProteinSolubility />
                </div>
            </div>
        </div>
    </div>
</template>
