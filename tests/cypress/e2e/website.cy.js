describe('Main page', () => {
    const inputPdbId = '#proteinFileInput';
    const inputAminoAcidSequenceId = '#inputSequence';
    const submitAminoAcidSequenceId = '#submitInference';
    const localisationImageClass = '.localisation-image';
    const foldingSelector = '#viewport div';

    it('Opened main page', () => {
        cy.visit('http://127.0.0.1:5000');
    });

    it('Amino Acid lab opened', () => {
        cy.visit('http://127.0.0.1:5000/amino-acid-page');
    });

    it('Amino-acid sequence input is visible', () => {
        cy.visit('http://127.0.0.1:5000/amino-acid-page');
        cy.get(inputAminoAcidSequenceId).should('exist');
    });

    it('Inference result is visible', () => {
        cy.visit('http://127.0.0.1:5000/amino-acid-page');
        cy.get(inputAminoAcidSequenceId).click();
        cy.get(inputAminoAcidSequenceId).type(`AAACGAGGCAA`);
        cy.get(submitAminoAcidSequenceId).click();
        cy.get(localisationImageClass).should('exist');
        cy.get(foldingSelector, {timeout: 90000}).should('exist');
    });

    //it('Drug target discovery lab opened', () => {
    //    cy.visit('http://127.0.0.1:5000/drug-target-discovery');
    //});
//
    //it('Pdb input is visible', () => {
    //    cy.visit('http://127.0.0.1:5000');
    //});
//
    //it('Inference result is visible', () => {
    //    cy.get('input[id=proteinFileInput]').selectFile('cypress/fixtures/test.pdb');
    //    cy.get('input[id=smilesFileInput]').selectFile('cypress/fixtures/test.sdf');
    //    cy.get('#submitInference').click();
    //    cy.get('#drugTargetDiscovery', {timeout: 90000}).should('exist');
    //});
})