<script lang="ts">
	export let name: string = "Roman";

	let rand = -1;
	function getRand() {
		fetch("./rand")
			.then((r) => r.text())
			.then((d) => (rand = d));
	}

	let file1;
	let file2;
	let files;
	$: if (file1 || file2) {
		// Note that `files` is of type `FileList`, not an Array:
		// https://developer.mozilla.org/en-US/docs/Web/API/FileList
		console.log(file1);
		console.log(file2);

		for (const file of file1) {
			console.log(`${file.name}: ${file.size} bytes`);
		}
		for (const file of file2) {
			console.log(`${file.name}: ${file.size} bytes`);
		}
	}

	async function doPost(e) {
		e.preventDefault();
		//const formData = new FormData(e.target);
		//formData.append("files", files);

		//const data = {};
		//for (let field of formData) {
		//	const [key, value] = field;
		//	data[key] = value;
		//}
		//console.log(data);

		const res = await fetch("./post", {
			method: "POST",
			//body: formData,
		});
	}
</script>

<main>
	<div>
		<h1>Hello {name}!</h1>
	</div>

	<div>
		<p>
			Your number is {rand}!
			<button on:click={getRand}>Get a random number</button>
		</p>
	</div>
	<div id="upload">
		<form on:submit|preventDefault={doPost}>
			<label for="upload_lower"> Upload a structure: </label>
			<input
				accept=".xyz, .cif, .in"
				bind:file1
				id="upload_lower"
				name="lower"
				type="file"
			/>
			<input
				accept=".xyz, .cif, .in"
				bind:file2
				id="upload_upper"
				name="upper"
				type="file"
			/>
			<button type="submit"> Upload </button>
		</form>
	</div>
</main>

{#if file1}
	<h2>Selected files:</h2>
	<p>{file1.length}</p>
	{#each Array.from(files) as file}
		<p>{file.name} ({file.size} bytes)</p>
	{/each}
{/if}

{#if file2}
	<h2>Selected files:</h2>
	<p>{file2.length}</p>
	{#each Array.from(files) as file}
		<p>{file.name} ({file.size} bytes)</p>
	{/each}
{/if}

<style>
	main {
		text-align: center;
		padding: 1em;
		max-width: 240px;
		margin: 0 auto;
	}

	h1 {
		color: #ff3e00;
		text-transform: uppercase;
		font-size: 4em;
		font-weight: 100;
	}

	* {
		box-sizing: border-box;
	}
	form {
		display: flex;
		flex-direction: column;
		width: 300px;
	}

	form > div {
		display: flex;
		justify-content: space-between;
	}

	form > div + * {
		margin-top: 10px;
	}

	@media (min-width: 640px) {
		main {
			max-width: none;
		}
	}
</style>
