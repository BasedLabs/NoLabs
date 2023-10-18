<script setup>
import store from '../../storage.js';
import { reactive } from 'vue';


const localisation = store.state.aminoAcid.inference.localisation;


const formatText = (prob) => {
    return prob.toFixed(2) * 100;
}

const highlightData = {
    mithochondria: {
        image: 'localisationImages/mithochondria.png',
        controlElementId: 'mithochondria-list-item',
        text: `Mithochondria ${formatText(localisation.mithochondria)}%`
    },
    nucleus: {
        image: 'localisationImages/nucleus.png',
        controlElementId: 'nucleus-list-item',
        text: `Nucleus ${formatText(localisation.nucleus)}%`
    },
    cytoplasm: {
        image: 'localisationImages/cytoplasm.png',
        controlElementId: 'cytoplasm-list-item',
        text: `Cytoplasm ${formatText(localisation.cytoplasm)}%`
    },
    other: {
        image: 'localisationImages/original.png',
        controlElementId: 'other-proteins-item',
        text: `Other proteins ${formatText(localisation.other)}%`
    },
    extracellular: {
        image: 'localisationImages/original.png',
        controlElementId: 'extracellular-proteins-item',
        text: `Extracellular ${formatText(localisation.extracellular)}%`
    },
    original: {
        image: 'localisationImages/original.png'
    }
};

const activePlace = reactive({
    activeObject: highlightData.original
})

const onMouseOver = (el) => {
    activePlace.activeObject = highlightData.original;

    let anyActive = false;
    for (const key in highlightData) {
        const data = highlightData[key];
        // If we are over the list item, render corresponding image
        if (data.controlElementId == el.target.id) {
            activePlace.activeObject = data;
            anyActive = true;
        }
    }

    // If not any active - render default
    if (!anyActive) {
        activePlace.activeObject = highlightData.original;
    }
}

const onMouseLeave = (el) => {
    activePlace.activeObject = highlightData.original;
}
</script>

<template>
    <div class="col-md-12" id="localisation">
        <img id="img" :src="activePlace.activeObject.image" class="localisation-image">
        <p>Hover on list</p>
        <ul class="list-group localisation-probs ">
            <li id="mithochondria-list-item" class="list-group-item"
                :class="highlightData.mithochondria === activePlace.activeObject ? 'active' : ''"
                v-on:mouseover="onMouseOver" v-on:mouseleave="onMouseLeave">{{ highlightData.mithochondria.text }}</li>
            <li id="nucleus-list-item" class="list-group-item"
                :class="highlightData.nucleus === activePlace.activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.nucleus.text }}</li>
            <li id="cytoplasm-list-item" class="list-group-item"
                :class="highlightData.cytoplasm === activePlace.activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.cytoplasm.text }}</li>
            <li id="other-proteins-item" class="list-group-item"
                :class="highlightData.other === activePlace.activeObject ? 'active' : ''" v-on:mouseover="onMouseOver"
                v-on:mouseleave="onMouseLeave">{{ highlightData.other.text }}</li>
            <li id="extracellular-proteins-item" class="list-group-item"
                :class="highlightData.extracellular === activePlace.activeObject ? 'active' : ''"
                v-on:mouseover="onMouseOver" v-on:mouseleave="onMouseLeave">{{ highlightData.extracellular.text }}</li>
        </ul>
    </div>
</template>
